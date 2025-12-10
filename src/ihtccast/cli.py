import argparse, os, numpy as np
from .core import SteelProps, SandProps, ModelConfig, h_parametric, run_identification

def main():
    p = argparse.ArgumentParser(description="IHTC identification for steel sand casting")
    p.add_argument("--temps", required=True, help="CSV with time_s and sensor columns")
    p.add_argument("--L_c", type=float, required=True, help="Casting thickness (m)")
    p.add_argument("--L_m", type=float, required=True, help="Mold thickness (m)")
    p.add_argument("--dx_c", type=float, default=5e-4, help="Casting dx (m)")
    p.add_argument("--dx_m", type=float, default=1e-3, help="Mold dx (m)")
    p.add_argument("--bc_mold", choices=["adiabatic","convection"], default="adiabatic")
    p.add_argument("--h_env", type=float, default=15.0)
    p.add_argument("--T_inf", type=float, default=298.0)
    p.add_argument("--t_end", type=float, required=True)
    p.add_argument("--dt_init", type=float, default=0.01)
    p.add_argument("--init", nargs="+", help="h_inf=..., h_peak=..., tau=..., t0=...")
    p.add_argument("--bounds", nargs="+", help="key=lo:hi ...")
    p.add_argument("--lambda", dest="lam", type=float, default=1e-3)
    p.add_argument("--bootstrap", type=int, default=0)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    os.makedirs(args.out, exist_ok=True)

    init_map = {kv.split("=")[0]: float(kv.split("=")[1]) for kv in (args.init or [])}
    bounds_map = {}
    for kv in (args.bounds or []):
        k, rng = kv.split("="); lo, hi = rng.split(":"); bounds_map[k]=(float(lo), float(hi))

    theta0 = [
        init_map.get('h_inf',800.0),
        init_map.get('h_peak',5000.0),
        init_map.get('tau',2.0),
        init_map.get('t0',0.0)
    ]
    lb = [
        bounds_map.get('h_inf',(50,5000))[0],
        bounds_map.get('h_peak',(200,20000))[0],
        bounds_map.get('tau',(0.02,30))[0],
        bounds_map.get('t0',(0.0,1.0))[0]
    ]
    ub = [
        bounds_map.get('h_inf',(50,5000))[1],
        bounds_map.get('h_peak',(200,20000))[1],
        bounds_map.get('tau',(0.02,30))[1],
        bounds_map.get('t0',(0.0,1.0))[1]
    ]

    cfg = ModelConfig(Lc=args.L_c, Lm=args.L_m, dx_c=args.dx_c, dx_m=args.dx_m,
                      T_init_cast=1750.0, T_init_mold=298.0,
                      bc_mold=args.bc_mold, h_env=args.h_env, T_inf=args.T_inf)

    results = run_identification(
        temps_csv=args.temps, cfg=cfg, theta0=np.array(theta0),
        bounds=(np.array(lb), np.array(ub)), lam=args.lam,
        n_boot=args.bootstrap, t_end=args.t_end, dt_init=args.dt_init, out_dir=args.out
    )
    print("Estimated parameters:", results["theta_hat"])
