import re

with open('crates/gneiss-rtk/src/engine/ppp_fg.rs', 'r') as f:
    content = f.read()

test_code = """
    #[test]
    fn test_solve_outlier_rejection_limit() {
        let fg = PppFactorGraph::default();
        let mut state = dummy_rtk_state();
        let mut sats = Vec::new();
        state.covariance = DMatrix::identity(CORE_STATE_SIZE + 5, CORE_STATE_SIZE + 5);
        for i in 1..=5 {
            let sat_id = SatelliteId { constellation: Constellation::Gps, prn: i };
            state.ambiguity_keys.push((sat_id, 0));
            state.ambiguities.push(0.0);
            
            // We need `sat_obs` to exist... but it's a reference! This makes creating loops hard.
        }
    }
"""

# Actually creating the sats in a loop with references is hard in Rust tests without an arena.
# But wait, what if I just use a test that doesn't need references?
