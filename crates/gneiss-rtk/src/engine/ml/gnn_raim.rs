use candle_core::{Result, Tensor, D};
use candle_nn::{linear, ops, Linear, Module, VarBuilder};

/// Multi-Layer Perceptron used for node embeddings and the output head
pub struct Mlp {
    fc1: Linear,
    fc2: Linear,
}

impl Mlp {
    pub fn new(in_dim: usize, hidden_dim: usize, out_dim: usize, vb: VarBuilder) -> Result<Self> {
        let fc1 = linear(in_dim, hidden_dim, vb.pp("fc1"))?;
        let fc2 = linear(hidden_dim, out_dim, vb.pp("fc2"))?;
        Ok(Self { fc1, fc2 })
    }

    pub fn forward(&self, xs: &Tensor) -> Result<Tensor> {
        let xs = self.fc1.forward(xs)?;
        let xs = xs.relu()?;
        self.fc2.forward(&xs)
    }
}

/// Simple Multi-Head Self-Attention to act as the GNN Message Passing layer
pub struct SelfAttention {
    q_proj: Linear,
    k_proj: Linear,
    v_proj: Linear,
    out_proj: Linear,
    embed_dim: usize,
    num_heads: usize,
}

impl SelfAttention {
    pub fn new(embed_dim: usize, num_heads: usize, vb: VarBuilder) -> Result<Self> {
        let q_proj = linear(embed_dim, embed_dim, vb.pp("q_proj"))?;
        let k_proj = linear(embed_dim, embed_dim, vb.pp("k_proj"))?;
        let v_proj = linear(embed_dim, embed_dim, vb.pp("v_proj"))?;
        let out_proj = linear(embed_dim, embed_dim, vb.pp("out_proj"))?;
        Ok(Self {
            q_proj,
            k_proj,
            v_proj,
            out_proj,
            embed_dim,
            num_heads,
        })
    }

    pub fn forward(&self, xs: &Tensor) -> Result<Tensor> {
        let (b_sz, seq_len) = (xs.dim(0)?, xs.dim(1)?);
        let (q, k, v) = self.project_qkv(xs, b_sz, seq_len)?;

        // Attention scores
        let att = q.matmul(&k.transpose(2, 3)?)?;
        let att = (att / ((self.embed_dim / self.num_heads) as f64).sqrt())?;
        let att = ops::softmax(&att, D::Minus1)?;

        let out = att
            .matmul(&v)?
            .transpose(1, 2)?
            .reshape((b_sz, seq_len, self.embed_dim))?;

        self.out_proj.forward(&out)
    }

    fn project_qkv(
        &self,
        xs: &Tensor,
        b_sz: usize,
        seq_len: usize,
    ) -> Result<(Tensor, Tensor, Tensor)> {
        let head_dim = self.embed_dim / self.num_heads;
        let q = self
            .q_proj
            .forward(xs)?
            .reshape((b_sz, seq_len, self.num_heads, head_dim))?
            .transpose(1, 2)?;
        let k = self
            .k_proj
            .forward(xs)?
            .reshape((b_sz, seq_len, self.num_heads, head_dim))?
            .transpose(1, 2)?;
        let v = self
            .v_proj
            .forward(xs)?
            .reshape((b_sz, seq_len, self.num_heads, head_dim))?
            .transpose(1, 2)?;
        Ok((q, k, v))
    }
}

/// The overall GNN RAIM Model
pub struct GnnRaimModel {
    embedding: Mlp,
    attention: SelfAttention,
    output_head: Mlp,
}

const FEATURE_DIM: usize = 5;
const HIDDEN_DIM: usize = 32;
const NUM_HEADS: usize = 4;
const HEAD_DIM: usize = HIDDEN_DIM / 2;
const OUT_DIM: usize = 1;

impl GnnRaimModel {
    /// Initializes the model structure (can be mocked with random weights via VarBuilder)
    pub fn new(vb: VarBuilder) -> Result<Self> {
        let embedding = Mlp::new(FEATURE_DIM, HIDDEN_DIM, HIDDEN_DIM, vb.pp("embedding"))?;
        // We use self-attention as the graph convolution mechanism
        let attention = SelfAttention::new(HIDDEN_DIM, NUM_HEADS, vb.pp("attention"))?;
        let output_head = Mlp::new(HIDDEN_DIM, HEAD_DIM, OUT_DIM, vb.pp("output_head"))?;

        Ok(Self {
            embedding,
            attention,
            output_head,
        })
    }

    /// Forward pass
    /// `node_features` expected shape: (batch_size, num_satellites, 5)
    pub fn forward(&self, node_features: &Tensor) -> Result<Tensor> {
        // 1. Independent node embedding
        let h = self.embedding.forward(node_features)?;

        // 2. Message passing (Self-Attention over the graph of satellites)
        let attn_out = self.attention.forward(&h)?;

        // Residual connection + non-linearity
        let h = (h + attn_out)?;
        let h = h.relu()?;

        // 3. Independent output projection for each node
        // logits shape: (batch_size, num_satellites, 1)
        let logits = self.output_head.forward(&h)?;

        // Output log-variance directly for numerical stability
        Ok(logits)
    }
}

/// Evaluates the GNN RAIM model on the current epoch's GNSS observations.
/// Returns a map of `SatelliteId` to predicted log-variance.
pub fn evaluate_gnn_raim(
    model: &GnnRaimModel,
    matched_obs: &[(crate::filter::DdObservation, crate::filter::DdObservation)],
    rov_llh: nalgebra::Vector3<f64>,
    pos_apc: nalgebra::Vector3<f64>,
    ephemerides: &[gneiss_core::ephemeris::Ephemeris],
    state_time: gneiss_core::time::GpsTime,
) -> std::collections::HashMap<gneiss_core::sat::SatelliteId, f64> {
    let mut map = std::collections::HashMap::new();
    let num_sats = matched_obs.len();
    if num_sats == 0 {
        return map;
    }

    let mut features = Vec::with_capacity(num_sats * FEATURE_DIM);
    let mut sats = Vec::with_capacity(num_sats);

    for (rov, _) in matched_obs {
        extract_and_push_features(
            rov,
            rov_llh,
            pos_apc,
            ephemerides,
            state_time,
            &mut features,
        );
        sats.push(rov.sat);
    }

    let dev = candle_core::Device::Cpu;
    if let Ok(tensor) = candle_core::Tensor::from_vec(features, (1, num_sats, FEATURE_DIM), &dev) {
        if let Ok(variances) = model.forward(&tensor) {
            if let Ok(vals) = variances.flatten_all() {
                if let Ok(vec) = vals.to_vec1::<f32>() {
                    for (i, v) in vec.into_iter().enumerate() {
                        let log_var = v as f64;
                        let variance = log_var.exp();
                        map.insert(sats[i], variance);
                    }
                }
            }
        }
    }

    map
}

/// Helper to extract features for a single observation and push to the flattened buffer
fn extract_and_push_features(
    rov: &crate::filter::DdObservation,
    rov_llh: nalgebra::Vector3<f64>,
    pos_apc: nalgebra::Vector3<f64>,
    ephemerides: &[gneiss_core::ephemeris::Ephemeris],
    state_time: gneiss_core::time::GpsTime,
    features: &mut Vec<f32>,
) {
    let mut el = 0.0;
    let mut az = 0.0;
    if let Some(eph) = ephemerides.iter().find(|e| e.sat() == rov.sat) {
        let (sat_pos, _) =
            crate::engine::measurement_math::get_sat_state(eph, rov.pr_l1, 0.0, state_time, pos_apc);
        let (a, e) = gneiss_core::coords::az_el(rov_llh, pos_apc, sat_pos);
        az = a.to_degrees();
        el = e.to_degrees();
    }
    let normalized = crate::engine::ml::dataset::normalize_features(
        rov.snr as f32,
        el as f32,
        az as f32,
        rov.doppler as f32,
        rov.locktime.unwrap_or(0) as f32,
    );

    features.extend_from_slice(&normalized);
}

#[cfg(test)]
mod tests {
    use super::*;
    use candle_core::{DType, Device, Tensor};
    use candle_nn::{VarBuilder, VarMap};

    #[test]
    fn test_gnn_raim_model_initialization_and_forward() -> Result<()> {
        let device = Device::Cpu;
        let varmap = VarMap::new();
        let vb = VarBuilder::from_varmap(&varmap, DType::F32, &device);
        let model = GnnRaimModel::new(vb)?;

        let batch_size = 1;
        let num_satellites = 10;

        let input = Tensor::randn(
            0f32,
            1f32,
            (batch_size, num_satellites, FEATURE_DIM),
            &device,
        )?;
        let output = model.forward(&input)?;

        assert_eq!(output.dims(), &[batch_size, num_satellites, OUT_DIM]);
        // The output is log-variance, so it can be any real number
        // We just ensure it doesn't return NaN here
        let vals = output.flatten_all()?.to_vec1::<f32>()?;
        for v in vals {
            assert!(!v.is_nan(), "Output log-variance must not be NaN");
        }

        Ok(())
    }
}
