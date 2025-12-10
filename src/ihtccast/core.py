import numpy as np, csv, os
from dataclasses import dataclass
from typing import Tuple, Dict
from scipy.optimize import least_squares
import matplotlib.pyplot as plt

@dataclass
class SteelProps:
    def k(self, T):  # W/m-K
        return 28.0 - 0.012*(T-300.0)
    def cp(self, T): # J/kg-K
        return 520.0 + 0.06*(T-300.0)
    def rho(self, T): # kg/m3
        return 7800.0*(1.0 - 1.2e-4*(T-300.0))
    Tsol: float = 1730.0  # K
    Tliq: float = 1800.0  # K
    L: float    = 2.8e5   # J/kg

@dataclass
class SandProps:
    def k(self, T):
        return 0.9 - 6e-4*(T-300.0)
    def cp(self, T):
        return 900.0 + 0.1*(T-300.0)
    def rho(self, T):
        return 1600.0

@dataclass
class Mesh1D:
    x: np.ndarray
    dx: float

@dataclass
class ModelConfig:
    Lc: float; Lm: float
    dx_c: float; dx_m: float
    T_init_cast: float; T_init_mold: float
    bc_mold: str = "adiabatic"
    h_env: float = 15.0
    T_inf: float = 298.0

def make_mesh(Lc, Lm, dx_c, dx_m):
    x_c = np.arange(0.0, Lc+1e-12, dx_c)
    x_m = np.arange(Lc+dx_m, Lc+Lm+1e-12, dx_m)
    return Mesh1D(x_c, dx_c), Mesh1D(x_m, dx_m)

def h_parametric(t, h_inf, h_peak, tau, t0):
    t = np.asarray(t)
    h = np.where(t>=t0, h_inf + (h_peak - h_inf)*np.exp(-(t - t0)/max(tau,1e-6)), h_inf)
    return np.clip(h, 1.0, None)

def read_temps_csv(path: str):
    with open(path, newline='') as f:
        reader = csv.DictReader(f)
        cols = reader.fieldnames
        rows = list(reader)
    arr = {k: np.array([float(r[k]) for r in rows]) for k in cols}
    t = arr['time_s']
    cast_cols = [c for c in cols if c.startswith("T_casting_")]
    mold_cols = [c for c in cols if c.startswith("T_mold_")]
    cast_depths = np.array([float(c.split("_")[2].replace("mm","")) for c in cast_cols])/1000.0
    mold_depths = np.array([float(c.split("_")[2].replace("mm","")) for c in mold_cols])/1000.0
    y_cast = np.vstack([arr[c] for c in cast_cols]).T if cast_cols else np.zeros((len(t),0))
    y_mold = np.vstack([arr[c] for c in mold_cols]).T if mold_cols else np.zeros((len(t),0))
    return t, cast_depths, mold_depths, y_cast, y_mold

def step_forward(cfg: ModelConfig, steel: SteelProps, sand: SandProps,
                 t_grid: np.ndarray, h_fun, T0_cast: float, T0_mold: float):
    mc, mm = make_mesh(cfg.Lc, cfg.Lm, cfg.dx_c, cfg.dx_m)
    Nc, Nm = len(mc.x), len(mm.x)
    Nt = len(t_grid)
    Tc = np.full(Nc, T0_cast, dtype=float)
    Tm = np.full(Nm, T0_mold, dtype=float)
    out_c = np.zeros((Nt, Nc)); out_m = np.zeros((Nt, Nm))
    out_c[0,:] = Tc; out_m[0,:] = Tm

    for n in range(1, Nt):
        dt = t_grid[n] - t_grid[n-1]
        kc = steel.k(Tc); rhoc = steel.rho(Tc); cpc = steel.cp(Tc)
        mush = (Tc>=steel.Tsol) & (Tc<=steel.Tliq)
        cpc_eff = cpc + steel.L*np.where(mush, 1.0/(steel.Tliq-steel.Tsol), 0.0)
        alpha_c = kc/(rhoc*cpc_eff)

        km = sand.k(Tm); rhom = sand.rho(Tm); cpm = sand.cp(Tm)
        alpha_m = km/(rhom*cpm)

        r_c = alpha_c*dt/(cfg.dx_c**2)
        r_m = alpha_m*dt/(cfg.dx_m**2)

        A_c = np.zeros((Nc, Nc)); b_c = Tc.copy()
        for i in range(Nc):
            if i==0:
                A_c[i,i]   = 1 + 2*r_c[i]
                A_c[i,i+1] = -2*r_c[i]
            elif i==Nc-1:
                A_c[i,i] = 1.0
            else:
                A_c[i,i-1] = -r_c[i]
                A_c[i,i]   = 1 + 2*r_c[i]
                A_c[i,i+1] = -r_c[i]

        A_m = np.zeros((Nm, Nm)); b_m = Tm.copy()
        for j in range(Nm):
            if j==Nm-1:
                if cfg.bc_mold=="adiabatic":
                    A_m[j,j]   = 1 + 2*r_m[j]
                    A_m[j,j-1] = -2*r_m[j]
                else:
                    beta = cfg.h_env*cfg.dx_m/km[j]
                    A_m[j,j]   = 1 + r_m[j]*(2+2*beta)
                    A_m[j,j-1] = -2*r_m[j]
                    b_m[j]    += 2*r_m[j]*beta*cfg.T_inf
            elif j==0:
                A_m[j,j] = 1.0
            else:
                A_m[j,j-1] = -r_m[j]
                A_m[j,j]   = 1 + 2*r_m[j]
                A_m[j,j+1] = -r_m[j]

        h = h_fun(t_grid[n])
        A_c[-1, -2] = -2*r_c[-1]
        A_c[-1, -1] = 1 + 2*r_c[-1]*(1 + h*cfg.dx_c/kc[-1])

        A_m[0,0] = 1 + 2*r_m[0]*(1 + h*cfg.dx_m/km[0])
        A_m[0,1] = -2*r_m[0]

        Tc_new = Tc.copy(); Tm_new = Tm.copy()
        for _ in range(3):
            rhs_c = b_c.copy()
            rhs_c[-1] += 2*r_c[-1]*h*cfg.dx_c/kc[-1]*Tm_new[0]
            Tc_new = np.linalg.solve(A_c, rhs_c)
            rhs_m = b_m.copy()
            rhs_m[0] += 2*r_m[0]*h*cfg.dx_m/km[0]*Tc_new[-1]
            Tm_new = np.linalg.solve(A_m, rhs_m)

        Tc, Tm = Tc_new, Tm_new
        out_c[n,:] = Tc; out_m[n,:] = Tm
    return out_c, out_m

def simulate_sensors(cfg, steel, sand, t, h_fun, cast_depths, mold_depths,
                     T0_cast, T0_mold):
    out_c, out_m = step_forward(cfg, steel, sand, t, h_fun, T0_cast, T0_mold)
    x_c = np.arange(0.0, cfg.Lc+1e-12, cfg.dx_c)
    x_m = np.arange(cfg.Lc+cfg.dx_m, cfg.Lc+cfg.Lm+1e-12, cfg.dx_m)
    y_cast = np.zeros((len(t), len(cast_depths)))
    y_mold = np.zeros((len(t), len(mold_depths)))
    for n in range(len(t)):
        for j,d in enumerate(cast_depths):
            xq = cfg.Lc - d
            idx = np.argmin(np.abs(x_c - xq))
            y_cast[n,j] = out_c[n, idx]
        for j,d in enumerate(mold_depths):
            xq = cfg.Lc + d
            idx = np.argmin(np.abs(x_m - xq))
            y_mold[n,j] = out_m[n, idx]
    return y_cast, y_mold

def residual_parametric(theta_vec, pack):
    h_inf, h_peak, tau, t0 = theta_vec
    def hfun(t): return h_parametric(t, h_inf, h_peak, tau, t0)
    y_cast_sim, y_mold_sim = simulate_sensors(**pack, h_fun=hfun)
    y_cast_meas, y_mold_meas = pack['y_cast_meas'], pack['y_mold_meas']
    w_cast, w_mold = pack['w_cast'], pack['w_mold']
    r_cast = (y_cast_sim - y_cast_meas).ravel()*w_cast
    r_mold = (y_mold_sim - y_mold_meas).ravel()*w_mold
    t = pack['t']
    h_vals = hfun(t)
    dh = np.diff(h_vals, prepend=h_vals[0])
    reg = pack['lam']*dh
    return np.concatenate([r_cast, r_mold, reg])

def run_identification(temps_csv, cfg: ModelConfig, theta0, bounds, lam, n_boot, t_end, dt_init, out_dir):
    steel = SteelProps(); sand = SandProps()
    t, cast_d, mold_d, y_cast_meas_raw, y_mold_meas_raw = read_temps_csv(temps_csv)
    t_fine = np.unique(np.concatenate([t, np.linspace(0, min(5.0,t_end), 200), np.linspace(5.0, t_end, 200)]))

    def interp_cols(Y): 
        return np.vstack([np.interp(t_fine, t, Y[:,j]) for j in range(Y.shape[1])]).T if Y.size else np.zeros((len(t_fine),0))
    y_cast_meas = interp_cols(y_cast_meas_raw)
    y_mold_meas = interp_cols(y_mold_meas_raw)

    w_vec = 1.0/np.sqrt(t_fine + 0.5)
    w_cast = np.repeat(w_vec[:,None], y_cast_meas.shape[1], axis=1).ravel() if y_cast_meas.size else np.array([])
    w_mold = np.repeat(w_vec[:,None], y_mold_meas.shape[1], axis=1).ravel() if y_mold_meas.size else np.array([])

    pack = dict(cfg=cfg, steel=steel, sand=sand, t=t_fine,
                cast_depths=cast_d, mold_depths=mold_d,
                T0_cast=cfg.T_init_cast, T0_mold=cfg.T_init_mold,
                y_cast_meas=y_cast_meas, y_mold_meas=y_mold_meas,
                w_cast=w_cast, w_mold=w_mold, lam=lam)

    res = least_squares(residual_parametric, theta0, args=(pack,),
                        bounds=bounds, xtol=1e-8, ftol=1e-8, gtol=1e-8, verbose=0)
    theta_hat = res.x

    def hfun(t): return h_parametric(t, *theta_hat)
    h_vals = hfun(t_fine)
    y_cast_sim, y_mold_sim = simulate_sensors(**pack, h_fun=hfun)

    os.makedirs(out_dir, exist_ok=True)
    np.savetxt(os.path.join(out_dir,"ihtc_curve.csv"), np.column_stack([t_fine, h_vals]),
               delimiter=",", header="t_s,h_Wm2K", comments="")

    plt.figure()
    plt.plot(t_fine, h_vals, label="IHTC estimate")
    plt.xlabel("Time (s)"); plt.ylabel("h (W m$^{-2}$ K$^{-1}$)")
    plt.legend(); plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "fig_h_curve.png"), dpi=200); plt.close()

    if y_cast_meas.size or y_mold_meas.size:
        plt.figure()
        if y_cast_meas.size:
            plt.plot(t_fine, y_cast_meas[:,0], label="Cast meas 1")
            plt.plot(t_fine, y_cast_sim[:,0], '--', label="Cast sim 1")
        if y_mold_meas.size:
            plt.plot(t_fine, y_mold_meas[:,0], label="Mold meas 1")
            plt.plot(t_fine, y_mold_sim[:,0], '--', label="Mold sim 1")
        plt.xlabel("Time (s)"); plt.ylabel("Temperature (K)")
        plt.legend(); plt.tight_layout()
        plt.savefig(os.path.join(out_dir,"fig_temp_compare.png"), dpi=200); plt.close()

    return {"theta_hat": theta_hat, "t": t_fine, "h": h_vals}
