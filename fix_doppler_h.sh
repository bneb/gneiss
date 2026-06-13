sed -i '' 's/h_row\[5\] = -e_los.z;/h_row\[5\] = -e_los.z;\n        if state_size > 16 {\n            h_row\[16\] = 1.0;\n        }/' crates/gneiss-rtk/src/doppler.rs
