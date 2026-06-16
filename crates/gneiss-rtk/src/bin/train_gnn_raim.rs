use candle_core::{Device, Result};
use candle_nn::{AdamW, Optimizer, ParamsAdamW, VarBuilder, VarMap};
use gneiss_rtk::engine::ml::gnn_raim::GnnRaimModel;
use gneiss_rtk::engine::ml::dataset_loader::DatasetLoader;
use gneiss_rtk::engine::ml::dataset::nll_loss_logvar;

fn main() -> Result<()> {
    let dataset_path = "shinjuku_gnn_dataset.csv";
    let device = Device::Cpu;

    println!("Loading dataset from {}...", dataset_path);
    let loader = DatasetLoader::new(dataset_path, &device)?;
    println!("Loaded {} epochs/batches.", loader.batches.len());

    let varmap = VarMap::new();
    let vb = VarBuilder::from_varmap(&varmap, candle_core::DType::F32, &device);
    let model = GnnRaimModel::new(vb)?;

    let params = ParamsAdamW {
        lr: 1e-4,
        weight_decay: 0.01,
        ..Default::default()
    };
    let mut optimizer = AdamW::new(varmap.all_vars(), params)?;

    let epochs = 10;
    use rand::seq::SliceRandom;
    let mut rng = rand::rng();

    println!("Starting training...");
    for epoch in 0..epochs {
        let mut epoch_loss = 0.0;
        let mut batches = 0;

        let mut shuffled_batches = loader.batches.clone();
        shuffled_batches.shuffle(&mut rng);

        for batch in &shuffled_batches {
            let features = &batch.features;
            let residuals = &batch.residuals;
            let mask = &batch.mask;

            let log_var = model.forward(features)?.squeeze(2)?;

            let loss = nll_loss_logvar(&log_var, residuals, mask)?;

            optimizer.backward_step(&loss)?;

            epoch_loss += loss.to_vec0::<f32>()?;
            batches += 1;
        }

        println!("Epoch {}: Mean Loss = {:.4}", epoch + 1, epoch_loss / (batches as f32));
    }

    varmap.save("gnn_raim_model.safetensors")?;
    println!("Model weights saved successfully to gnn_raim_model.safetensors");

    Ok(())
}
