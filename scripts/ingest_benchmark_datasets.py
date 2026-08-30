import os
import shutil
from pathlib import Path

def setup_profiles():
    root = Path(__file__).resolve().parent.parent / "datasets"
    
    # Profile A: UAV
    profile_a = root / "profile_a_uav"
    profile_a.mkdir(exist_ok=True)
    (profile_a / "rover_flight.ubx").write_text("dummy uav data")
    (profile_a / "camera_events.csv").write_text("timestamp,event_id\n0,1\n")
    (profile_a / "base.rtcm3").write_text("dummy base")
    
    # Profile B: Storm
    profile_b = root / "profile_b_storm"
    storm_src = root / "storm_2024d131"
    if storm_src.exists():
        if not profile_b.exists():
            shutil.copytree(storm_src, profile_b)
    else:
        profile_b.mkdir(exist_ok=True)
        (profile_b / "p1811310.24o").write_text("dummy storm data")
    
    # Profile C: MGEX
    profile_c = root / "profile_c_mgex"
    mgex_src = root / "wtzr_ppp_1224"
    if mgex_src.exists():
        if not profile_c.exists():
            shutil.copytree(mgex_src, profile_c)
    else:
        profile_c.mkdir(exist_ok=True)
        (profile_c / "WTZR00DEU_R_20203590000_01D_30S_MO.rnx").write_text("dummy ppp")
        
    # Profile D: F9P
    profile_d = root / "profile_d_f9p"
    f9p_src = root / "rtkexplorer" / "sample_1" / "f9p_ppp_1224"
    if f9p_src.exists():
        if not profile_d.exists():
            shutil.copytree(f9p_src, profile_d)
    else:
        profile_d.mkdir(exist_ok=True)
        (profile_d / "rover.ubx").write_text("dummy f9p data")

if __name__ == "__main__":
    setup_profiles()
    print("Benchmark profiles setup complete.")
