sed -i '' 's/} else {/} else if j == 16 {\n                    assert!((h\[(i, j)\] - 1.0).abs() < 1e-6, "H[{},16] should be 1.0", i);\n                } else {/' crates/gneiss-rtk/src/doppler.rs
