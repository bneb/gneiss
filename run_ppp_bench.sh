#!/bin/bash
cargo build --release
target/release/gneiss-cli process --mode ppp --rover datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover.obs --nav datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover.nav --output out_ppp.pos --sp3 datasets/rtkexplorer/sample_1/f9p_ppp_1224/ESA0MGNFIN_20203590000_01D_05M_ORB.SP3 --clk datasets/rtkexplorer/sample_1/f9p_ppp_1224/ESA0MGNFIN_20203590000_01D_30S_CLK.CLK --config datasets/rtkexplorer/sample_1/f9p_ppp_1224/f9p_config.json
target/release/gneiss-cli eval --solution out_ppp.pos --truth datasets/rtkexplorer/sample_1/f9p_ppp_1224/rover_ppk.pos
