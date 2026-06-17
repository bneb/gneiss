with open("crates/gneiss-rtk/src/engine/fgo/solver.rs", "r") as f:
    content = f.read()

content = content.replace(
"""        let success = lm.optimize(&mut graph);
        assert!(success); // it should hit max_iterations and fail to converge because residual doesn't shrink.
    }""",
"""        let success = lm.optimize(&mut graph);
        assert!(success);
        
        if let Some(VariableType::Ambiguity(a)) = graph.variables.get(&1) {
            assert!((a - 5.0).abs() < 1e-3, "State did not converge, got {}", a);
        } else {
            panic!("Variable not found");
        }
    }""")

with open("crates/gneiss-rtk/src/engine/fgo/solver.rs", "w") as f:
    f.write(content)
