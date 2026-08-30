#![allow(clippy::unwrap_used)]

use super::*;


    #[test]
    fn test_nhc_factor_forward_motion_zero_lateral_residual() {
        let v_pose = VariableId::new(1);
        let v_vel = VariableId::new(2);

        let mut graph_vars = std::collections::BTreeMap::new();
        graph_vars.insert(v_pose, crate::swfg::variables::VariableNode {
            id: v_pose,
            kind: crate::swfg::variables::VariableKind::Pose { epoch: 0 },
            value: nalgebra::DVector::from_vec(vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]),
        });
        graph_vars.insert(v_vel, crate::swfg::variables::VariableNode {
            id: v_vel,
            kind: crate::swfg::variables::VariableKind::Velocity { epoch: 0 },
            value: nalgebra::DVector::from_vec(vec![10.0, 0.0, 0.0]),
        });

        let values = VariableValues::build(&graph_vars);

        let nhc = NhcFactor {
            var_pose: v_pose,
            var_vel: v_vel,
            variance_y: 0.01,
            variance_z: 0.01,
            variables: vec![v_pose, v_vel],
        };

        let r = nhc.residual(&values);
        assert_eq!(r.len(), 2);
        assert!(r[0].abs() < 1e-6, "lateral velocity should be 0, got {}", r[0]);
        assert!(r[1].abs() < 1e-6, "vertical velocity should be 0, got {}", r[1]);
    }
