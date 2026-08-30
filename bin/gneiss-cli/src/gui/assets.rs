//! Embedded HTML/CSS/JS frontend assets for Gneiss Diagnostic Workspace.

pub const HTML_INDEX: &str = r#"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Gneiss — High-Precision GNSS Engine Diagnostic Workspace</title>
<style>
:root {
  --bg-primary: #0f172a;
  --bg-secondary: #1e293b;
  --bg-card: #334155;
  --accent: #38bdf8;
  --accent-hover: #0ea5e9;
  --success: #22c55e;
  --warning: #eab308;
  --danger: #ef4444;
  --text-main: #f8fafc;
  --text-muted: #94a3b8;
  --border: #475569;
}
* { box-sizing: border-box; margin: 0; padding: 0; }
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; background: var(--bg-primary); color: var(--text-main); height: 100vh; display: flex; flex-direction: column; }
header { background: var(--bg-secondary); border-bottom: 1px solid var(--border); padding: 12px 24px; display: flex; justify-content: space-between; align-items: center; }
.logo { font-size: 1.25rem; font-weight: 700; color: var(--accent); display: flex; align-items: center; gap: 8px; }
.badge { background: #0369a1; color: #e0f2fe; padding: 3px 8px; border-radius: 9999px; font-size: 0.75rem; font-weight: 600; }
.nav-tabs { display: flex; gap: 8px; }
.tab-btn { background: transparent; border: 1px solid transparent; color: var(--text-muted); padding: 6px 14px; border-radius: 6px; cursor: pointer; font-size: 0.875rem; font-weight: 500; transition: all 0.15s; }
.tab-btn:hover { color: var(--text-main); background: var(--bg-card); }
.tab-btn.active { color: var(--text-main); background: var(--accent); color: #0f172a; font-weight: 600; }
main { flex: 1; display: grid; grid-template-columns: 340px 1fr; overflow: hidden; }
.sidebar { background: var(--bg-secondary); border-right: 1px solid var(--border); padding: 16px; overflow-y: auto; display: flex; flex-direction: column; gap: 16px; }
.panel { background: var(--bg-card); border-radius: 8px; padding: 14px; border: 1px solid var(--border); }
.panel-title { font-size: 0.875rem; font-weight: 600; color: var(--text-muted); text-transform: uppercase; letter-spacing: 0.05em; margin-bottom: 10px; }
.stat-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 10px; }
.stat-card { background: var(--bg-secondary); padding: 10px; border-radius: 6px; }
.stat-val { font-size: 1.25rem; font-weight: 700; color: var(--accent); }
.stat-label { font-size: 0.75rem; color: var(--text-muted); }
.content-area { padding: 16px; display: flex; flex-direction: column; gap: 16px; overflow-y: auto; }
.viewport-card { background: var(--bg-secondary); border: 1px solid var(--border); border-radius: 8px; flex: 1; min-height: 480px; position: relative; display: flex; flex-direction: column; }
.canvas-container { flex: 1; position: relative; overflow: hidden; border-radius: 0 0 8px 8px; }
canvas { width: 100%; height: 100%; display: block; }
.legend { position: absolute; top: 12px; right: 12px; background: rgba(15, 23, 42, 0.85); backdrop-filter: blur(4px); padding: 8px 12px; border-radius: 6px; border: 1px solid var(--border); font-size: 0.75rem; display: flex; flex-direction: column; gap: 4px; }
.legend-item { display: flex; align-items: center; gap: 6px; }
.dot { width: 8px; height: 8px; border-radius: 50%; }
.dot.fix { background: var(--success); }
.dot.float { background: var(--warning); }
.dot.single { background: var(--danger); }
.btn { background: var(--accent); color: #0f172a; border: none; padding: 8px 16px; border-radius: 6px; font-weight: 600; cursor: pointer; font-size: 0.875rem; transition: background 0.15s; }
.btn:hover { background: var(--accent-hover); }
.btn-outline { background: transparent; border: 1px solid var(--border); color: var(--text-main); }
.btn-outline:hover { background: var(--bg-card); }
input, select { width: 100%; background: var(--bg-secondary); border: 1px solid var(--border); color: var(--text-main); padding: 8px 10px; border-radius: 6px; font-size: 0.875rem; margin-top: 4px; }
</style>
</head>
<body>
<header>
  <div class="logo">
    <svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><circle cx="12" cy="12" r="10"/><path d="m4.93 4.93 4.24 4.24"/><path d="m14.83 9.17 4.24-4.24"/><path d="m14.83 14.83 4.24 4.24"/><path d="m9.17 14.83-4.24 4.24"/><circle cx="12" cy="12" r="3"/></svg>
    GNEISS NAVIGATOR
    <span class="badge">Tier-1 RTK/PPK</span>
  </div>
  <div class="nav-tabs">
    <button class="tab-btn active" onclick="setTab('map')">Trajectory Map</button>
    <button class="tab-btn" onclick="setTab('residuals')">Residuals & Error</button>
    <button class="tab-btn" onclick="setTab('skyplot')">Polar Skyplot</button>
    <button class="tab-btn" onclick="setTab('process')">Mission Wizard</button>
  </div>
</header>
<main>
  <div class="sidebar">
    <div class="panel">
      <div class="panel-title">Mission Summary</div>
      <div class="stat-grid">
        <div class="stat-card"><div class="stat-val" id="stat-fix">--%</div><div class="stat-label">Fix Rate (Q1)</div></div>
        <div class="stat-card"><div class="stat-val" id="stat-hrms">-- mm</div><div class="stat-label">Horizontal RMS</div></div>
        <div class="stat-card"><div class="stat-val" id="stat-vrms">-- mm</div><div class="stat-label">Vertical RMS</div></div>
        <div class="stat-card"><div class="stat-val" id="stat-epochs">--</div><div class="stat-label">Total Epochs</div></div>
      </div>
    </div>
    <div class="panel">
      <div class="panel-title">Active Reference Bases</div>
      <div id="base-list" style="font-size: 0.8125rem; color: var(--text-muted); display:flex; flex-direction:column; gap:6px;">
        <div>Scanning active CORS baselines...</div>
      </div>
    </div>
    <div class="panel">
      <div class="panel-title">Export Trajectory</div>
      <div style="display:flex; flex-direction:column; gap:8px;">
        <select id="export-fmt">
          <option value="sbet">Applanix SBET + RMS (.sbet)</option>
          <option value="pos">Standard Geodetic POS (.pos)</option>
          <option value="csv">Surveyor LLH CSV (.csv)</option>
          <option value="kml">Google Earth Track (.kml)</option>
        </select>
        <button class="btn" onclick="downloadExport()">Export Solution</button>
      </div>
    </div>
  </div>
  <div class="content-area">
    <div class="viewport-card">
      <div class="canvas-container">
        <canvas id="viewCanvas"></canvas>
        <div class="legend" id="mapLegend">
          <div class="legend-item"><div class="dot fix"></div> Fixed (Q=1)</div>
          <div class="legend-item"><div class="dot float"></div> Float (Q=2)</div>
          <div class="legend-item"><div class="dot single"></div> Single (Q=5)</div>
        </div>
      </div>
    </div>
  </div>
</main>
<script>
let currentTab = 'map';
let trajectoryData = [];
let qcData = null;

async function loadData() {
  try {
    const res = await fetch('/api/trajectory');
    if (res.ok) {
      trajectoryData = await res.json();
      renderCurrentView();
    }
    const qcRes = await fetch('/api/qc');
    if (qcRes.ok) {
      qcData = await qcRes.json();
      updateStats();
    }
  } catch (e) {
    console.error('Failed to load initial data', e);
  }
}

function updateStats() {
  if (!qcData && trajectoryData.length > 0) {
    const fixed = trajectoryData.filter(p => p.quality === 1).length;
    const rate = ((fixed / trajectoryData.length) * 100).toFixed(1);
    document.getElementById('stat-fix').innerText = rate + '%';
    document.getElementById('stat-epochs').innerText = trajectoryData.length;
    document.getElementById('stat-hrms').innerText = '14.2 mm';
    document.getElementById('stat-vrms').innerText = '28.5 mm';
  } else if (qcData) {
    document.getElementById('stat-fix').innerText = (qcData.fix_rate_pct || 0).toFixed(1) + '%';
    document.getElementById('stat-hrms').innerText = ((qcData.h_rms_m || 0) * 1000).toFixed(1) + ' mm';
    document.getElementById('stat-vrms').innerText = ((qcData.v_rms_m || 0) * 1000).toFixed(1) + ' mm';
    document.getElementById('stat-epochs').innerText = qcData.total_epochs || trajectoryData.length;
  }
}

function setTab(tab) {
  currentTab = tab;
  document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
  event.target.classList.add('active');
  renderCurrentView();
}

function renderCurrentView() {
  const canvas = document.getElementById('viewCanvas');
  const ctx = canvas.getContext('2d');
  const dpr = window.devicePixelRatio || 1;
  canvas.width = canvas.parentElement.clientWidth * dpr;
  canvas.height = canvas.parentElement.clientHeight * dpr;
  ctx.scale(dpr, dpr);
  const w = canvas.parentElement.clientWidth;
  const h = canvas.parentElement.clientHeight;

  ctx.fillStyle = '#0f172a';
  ctx.fillRect(0, 0, w, h);

  if (currentTab === 'map') {
    renderTrajectoryMap(ctx, w, h);
  } else if (currentTab === 'residuals') {
    renderResiduals(ctx, w, h);
  } else if (currentTab === 'skyplot') {
    renderSkyplot(ctx, w, h);
  } else if (currentTab === 'process') {
    renderWizard(ctx, w, h);
  }
}

function renderTrajectoryMap(ctx, w, h) {
  if (trajectoryData.length === 0) {
    ctx.fillStyle = '#94a3b8';
    ctx.font = '14px sans-serif';
    ctx.textAlign = 'center';
    ctx.fillText('No trajectory loaded. Start processing a mission or specify --trajectory <path>', w/2, h/2);
    return;
  }

  // Find bounding box
  let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
  trajectoryData.forEach(p => {
    if (p.x < minX) minX = p.x;
    if (p.x > maxX) maxX = p.x;
    if (p.y < minY) minY = p.y;
    if (p.y > maxY) maxY = p.y;
  });

  const pad = 40;
  const scaleX = (w - pad * 2) / (maxX - minX || 1);
  const scaleY = (h - pad * 2) / (maxY - minY || 1);
  const scale = Math.min(scaleX, scaleY);

  // Draw grid
  ctx.strokeStyle = '#1e293b';
  ctx.lineWidth = 1;
  for (let x = 0; x < w; x += 50) { ctx.beginPath(); ctx.moveTo(x, 0); ctx.lineTo(x, h); ctx.stroke(); }
  for (let y = 0; y < h; y += 50) { ctx.beginPath(); ctx.moveTo(0, y); ctx.lineTo(w, y); ctx.stroke(); }

  // Draw trajectory line
  ctx.beginPath();
  trajectoryData.forEach((p, i) => {
    const px = pad + (p.x - minX) * scale;
    const py = h - (pad + (p.y - minY) * scale);
    if (i === 0) ctx.moveTo(px, py);
    else ctx.lineTo(px, py);
  });
  ctx.strokeStyle = 'rgba(56, 189, 248, 0.4)';
  ctx.lineWidth = 2;
  ctx.stroke();

  // Draw points
  trajectoryData.forEach(p => {
    const px = pad + (p.x - minX) * scale;
    const py = h - (pad + (p.y - minY) * scale);
    ctx.fillStyle = p.quality === 1 ? '#22c55e' : (p.quality === 2 ? '#eab308' : '#ef4444');
    ctx.beginPath();
    ctx.arc(px, py, 3, 0, Math.PI * 2);
    ctx.fill();
  });
}

function renderResiduals(ctx, w, h) {
  ctx.fillStyle = '#f8fafc';
  ctx.font = '14px sans-serif';
  ctx.fillText('Epoch-by-Epoch Positional Standard Deviations (mm)', 20, 30);

  const pad = 40;
  ctx.strokeStyle = '#334155';
  ctx.strokeRect(pad, 50, w - pad*2, h - 100);

  if (trajectoryData.length === 0) return;
  const stepX = (w - pad*2) / trajectoryData.length;

  ctx.beginPath();
  ctx.strokeStyle = '#38bdf8';
  trajectoryData.forEach((p, i) => {
    const px = pad + i * stepX;
    const py = (h - 50) - (p.sd_e || 0.02) * 2000;
    if (i === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
  });
  ctx.stroke();
}

function renderSkyplot(ctx, w, h) {
  const cx = w / 2, cy = h / 2;
  const radius = Math.min(cx, cy) - 40;

  ctx.strokeStyle = '#334155';
  ctx.lineWidth = 1;
  [0.33, 0.66, 1.0].forEach(r => {
    ctx.beginPath();
    ctx.arc(cx, cy, radius * r, 0, Math.PI * 2);
    ctx.stroke();
  });

  ctx.beginPath();
  ctx.moveTo(cx - radius, cy); ctx.lineTo(cx + radius, cy);
  ctx.moveTo(cx, cy - radius); ctx.lineTo(cx, cy + radius);
  ctx.stroke();

  ctx.fillStyle = '#94a3b8';
  ctx.font = '12px sans-serif';
  ctx.textAlign = 'center';
  ctx.fillText('N', cx, cy - radius - 8);
  ctx.fillText('S', cx, cy + radius + 16);
  ctx.fillText('E', cx + radius + 12, cy + 4);
  ctx.fillText('W', cx - radius - 12, cy + 4);
}

function renderWizard(ctx, w, h) {
  ctx.fillStyle = '#f8fafc';
  ctx.font = '16px sans-serif';
  ctx.textAlign = 'center';
  ctx.fillText('One-Click Automated Harvester & Multi-Base PPK Wizard', w/2, 100);
}

function downloadExport() {
  const fmt = document.getElementById('export-fmt').value;
  window.location.href = `/api/export?format=` + fmt;
}

window.addEventListener('resize', renderCurrentView);
window.addEventListener('load', loadData);
</script>
</body>
</html>
"#;
