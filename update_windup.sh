sed -i '' 's/let prev_w_bas_sat = \*state.windup.get(&ctx.base_sat.sat).unwrap_or(&0.0);/let prev_w_bas_sat = \*state.base_windup.get(\&ctx.base_sat.sat).unwrap_or(\&0.0);/' crates/gneiss-rtk/src/engine/measurement.rs
sed -i '' 's/let prev_w_bas_ref = \*state.windup.get(&ctx.ref_base.sat).unwrap_or(&0.0);/let prev_w_bas_ref = \*state.base_windup.get(\&ctx.ref_base.sat).unwrap_or(\&0.0);/' crates/gneiss-rtk/src/engine/measurement.rs
sed -i '' 's/state.windup.insert(ctx.base_sat.sat, w_bas_sat);/state.base_windup.insert(ctx.base_sat.sat, w_bas_sat);/' crates/gneiss-rtk/src/engine/measurement.rs
sed -i '' 's/state.windup.insert(ctx.ref_base.sat, w_bas_ref);/state.base_windup.insert(ctx.ref_base.sat, w_bas_ref);/' crates/gneiss-rtk/src/engine/measurement.rs
