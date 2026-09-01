//! Embedded HTML/CSS/JS frontend assets for Gneiss Diagnostic Workspace.

pub const HTML_INDEX: &str = r#"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"><meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Gneiss — High-Precision GNSS Engine Diagnostic Workspace</title>
<link rel="preconnect" href="https://fonts.googleapis.com"><link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Newsreader:ital,opsz,wght@0,6..72,400..700;1,6..72,400..700&family=JetBrains+Mono:wght@400;500;600;700&family=Inter:wght@400;500;600;700&display=swap" rel="stylesheet">
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css" crossorigin=""/>
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js" crossorigin=""></script>
<style>
:root {
  --bg-app: #080c14; --bg-surface: #101726; --bg-elevated: #182236;
  --accent: #38bdf8; --success: #10b981; --warning: #f59e0b; --danger: #ef4444;
  --text-pure: #ffffff; --text-main: #f1f5f9; --text-muted: #94a3b8;
  --border: #233148; --font-serif: 'Newsreader', Charter, Georgia, serif;
  --font-sans: 'Inter', -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
  --font-mono: 'JetBrains Mono', ui-monospace, SFMono-Regular, Menlo, monospace;
}
* { box-sizing: border-box; margin: 0; padding: 0; }
body { font-family: var(--font-sans); background: var(--bg-app); color: var(--text-main); height: 100vh; display: flex; flex-direction: column; overflow: hidden; -webkit-font-smoothing: antialiased; }
header { background: var(--bg-surface); border-bottom: 1px solid var(--border); padding: 10px 24px; display: flex; justify-content: space-between; align-items: center; z-index: 1000; box-shadow: 0 2px 8px rgba(0,0,0,0.4); }
.brand { display: flex; align-items: baseline; gap: 8px; }
.brand-title { font-family: var(--font-serif); font-size: 1.65rem; font-style: italic; font-weight: 500; color: var(--text-pure); letter-spacing: -0.01em; line-height: 1; }
.brand-sub { font-size: 0.7rem; font-weight: 700; color: var(--accent); letter-spacing: 0.15em; text-transform: uppercase; }
.badge { background: #0c2b45; color: #7dd3fc; border: 1px solid #0284c7; padding: 2px 8px; border-radius: 4px; font-size: 0.6875rem; font-weight: 600; font-family: var(--font-mono); }
.nav-tabs { display: flex; background: var(--bg-app); padding: 3px; border-radius: 6px; border: 1px solid var(--border); gap: 2px; }
.tab-btn { background: transparent; border: none; color: var(--text-muted); padding: 6px 14px; border-radius: 4px; cursor: pointer; font-size: 0.8125rem; font-weight: 600; transition: all 0.15s ease; }
.tab-btn:hover { color: var(--text-pure); }
.tab-btn.active { background: var(--bg-elevated); color: var(--accent); border: 1px solid var(--border); box-shadow: inset 0 1px 0 rgba(255,255,255,0.08); }
.wizard-btn { background: #0369a1; color: #fff; border: 1px solid var(--accent); padding: 6px 12px; border-radius: 4px; font-size: 0.75rem; font-weight: 700; cursor: pointer; }
.wizard-btn:hover { background: #0284c7; }
main { flex: 1; display: grid; grid-template-columns: 310px 1fr; overflow: hidden; }
.sidebar { background: var(--bg-surface); border-right: 1px solid var(--border); padding: 14px; overflow-y: auto; display: flex; flex-direction: column; gap: 12px; z-index: 999; }
.instrument-card { background: var(--bg-elevated); border: 1px solid var(--border); border-radius: 6px; padding: 12px; }
.card-header { font-size: 0.7rem; font-weight: 700; color: var(--text-muted); text-transform: uppercase; letter-spacing: 0.1em; margin-bottom: 8px; display: flex; justify-content: space-between; }
.telemetry-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 6px; }
.telemetry-box { background: var(--bg-app); border: 1px solid var(--border); padding: 8px; border-radius: 4px; }
.telemetry-val { font-family: var(--font-mono); font-size: 1.15rem; font-weight: 700; color: var(--accent); font-variant-numeric: tabular-nums; }
.telemetry-lbl { font-size: 0.6875rem; color: var(--text-muted); margin-top: 2px; }
.base-item { background: var(--bg-app); border: 1px solid var(--border); padding: 6px 8px; border-radius: 4px; font-family: var(--font-mono); font-size: 0.75rem; cursor: pointer; transition: border-color 0.15s; }
.base-item:hover { border-color: var(--accent); background: #131f33; }
.base-name { font-weight: 700; color: #38bdf8; display: flex; align-items: center; gap: 5px; }
.content-area { padding: 12px; display: flex; flex-direction: column; gap: 12px; overflow: hidden; position: relative; }
.viewport-frame { background: var(--bg-surface); border: 1px solid var(--border); border-radius: 6px; flex: 1; position: relative; display: flex; flex-direction: column; overflow: hidden; }
#mapContainer, #chartCanvas { position: absolute; top: 0; left: 0; width: 100%; height: 100%; background: var(--bg-app); }
.map-switch { position: absolute; top: 10px; left: 52px; z-index: 999; background: var(--bg-surface); border: 1px solid var(--border); border-radius: 6px; padding: 4px; display: flex; gap: 4px; box-shadow: 0 4px 12px rgba(0,0,0,0.5); }
.switch-btn { background: transparent; border: 1px solid transparent; color: var(--text-muted); padding: 4px 10px; border-radius: 4px; font-size: 0.75rem; font-weight: 600; cursor: pointer; }
.switch-btn.active { background: var(--bg-elevated); color: var(--accent); border-color: var(--border); }
.leaflet-bar { border: 1px solid var(--border) !important; box-shadow: 0 4px 12px rgba(0,0,0,0.5) !important; border-radius: 6px !important; overflow: hidden; }
.leaflet-bar a { background-color: var(--bg-surface) !important; color: var(--text-main) !important; border-bottom: 1px solid var(--border) !important; }
.leaflet-bar a:hover { background-color: var(--bg-elevated) !important; color: var(--accent) !important; }
.legend-box { position: absolute; top: 10px; right: 10px; z-index: 999; background: var(--bg-surface); border: 1px solid var(--border); border-radius: 6px; padding: 8px 12px; font-size: 0.75rem; display: flex; flex-direction: column; gap: 4px; font-family: var(--font-mono); }
.legend-row { display: flex; align-items: center; gap: 6px; }
.jewel { width: 8px; height: 8px; border-radius: 2px; }
.jewel.fix { background: var(--success); box-shadow: 0 0 6px var(--success); }
.jewel.float { background: var(--warning); }
.jewel.single { background: var(--danger); }
.jewel.base { background: var(--accent); }
.playback-bar { position: absolute; bottom: 12px; left: 12px; right: 12px; z-index: 999; background: var(--bg-surface); border: 1px solid var(--border); border-radius: 6px; padding: 8px 14px; display: flex; align-items: center; gap: 12px; font-family: var(--font-mono); }
.play-btn { background: var(--bg-elevated); border: 1px solid var(--border); color: var(--accent); width: 32px; height: 32px; border-radius: 4px; cursor: pointer; font-weight: 700; display: flex; align-items: center; justify-content: center; }
.scrub-slider { flex: 1; accent-color: var(--accent); cursor: pointer; height: 6px; }
.scrub-time { font-size: 0.75rem; color: var(--text-main); white-space: nowrap; }
.action-btn { background: #0284c7; color: #ffffff; border: 1px solid #38bdf8; padding: 8px 12px; border-radius: 4px; font-weight: 600; cursor: pointer; font-size: 0.8125rem; width: 100%; }
select, input[type="text"] { width: 100%; background: var(--bg-app); border: 1px solid var(--border); color: var(--text-main); padding: 7px 9px; border-radius: 4px; font-size: 0.75rem; font-family: var(--font-mono); }
.modal-backdrop { position: fixed; top: 0; left: 0; width: 100vw; height: 100vh; background: rgba(8,12,20,0.85); z-index: 2000; display: none; align-items: center; justify-content: center; }
.modal-box { background: var(--bg-surface); border: 1px solid var(--border); border-radius: 8px; width: 520px; max-width: 90vw; padding: 20px; display: flex; flex-direction: column; gap: 14px; }
.modal-header { display: flex; justify-content: space-between; align-items: center; border-bottom: 1px solid var(--border); padding-bottom: 10px; }
.modal-title { font-family: var(--font-serif); font-size: 1.35rem; font-style: italic; font-weight: 600; color: #fff; }
.modal-close { background: transparent; border: none; color: var(--text-muted); font-size: 1.25rem; cursor: pointer; }
.form-group { display: flex; flex-direction: column; gap: 4px; }
.form-label { font-size: 0.7rem; font-weight: 700; color: var(--text-muted); text-transform: uppercase; }
.err-box { background: #450a0a; border: 1px solid #dc2626; color: #fca5a5; padding: 8px; border-radius: 4px; font-size: 0.75rem; font-family: var(--font-mono); display: none; }
</style>
</head>
<body>
<header>
  <div class="brand"><span class="brand-title">Gneiss</span><span class="brand-sub">Navigator</span><span class="badge">TIER-1 PRECISION</span></div>
  <div class="nav-tabs">
    <button class="tab-btn active" onclick="setTab('map')">Geospatial Map</button>
    <button class="tab-btn" onclick="setTab('residuals')">Residuals & Error (mm)</button>
    <button class="tab-btn" onclick="setTab('skyplot')">Polar Skyplot</button>
  </div>
  <button class="wizard-btn" onclick="openWizard()">+ New Mission Wizard</button>
</header>
<main>
  <div class="sidebar">
    <div class="instrument-card">
      <div class="card-header"><span>Telemetry Summary</span><span style="color:var(--success)">ONLINE</span></div>
      <div class="telemetry-grid">
        <div class="telemetry-box"><div class="telemetry-val" id="stat-fix">--%</div><div class="telemetry-lbl">Fix Rate (Q1)</div></div>
        <div class="telemetry-box"><div class="telemetry-val" id="stat-hrms">-- mm</div><div class="telemetry-lbl">Horizontal RMS</div></div>
        <div class="telemetry-box"><div class="telemetry-val" id="stat-vrms">-- mm</div><div class="telemetry-lbl">Vertical RMS</div></div>
        <div class="telemetry-box"><div class="telemetry-val" id="stat-epochs">--</div><div class="telemetry-lbl">Total Epochs</div></div>
      </div>
    </div>
    <div class="instrument-card">
      <div class="card-header">Reference Base Stations</div>
      <div id="base-list" style="display:flex; flex-direction:column; gap:6px;"><div>Scanning active bases...</div></div>
    </div>
    <div class="instrument-card">
      <div class="card-header">Export Trajectory</div>
      <div style="display:flex; flex-direction:column; gap:8px;">
        <select id="export-fmt">
          <option value="sbet">Applanix POSPac SBET (.sbet)</option><option value="pos">Standard Geodetic POS (.pos)</option>
          <option value="csv">Surveyor LLH CSV (.csv)</option><option value="kml">Google Earth Track (.kml)</option>
        </select>
        <button class="action-btn" onclick="downloadExport()">Export Solution</button>
      </div>
    </div>
  </div>
  <div class="content-area">
    <div class="viewport-frame">
      <div id="mapContainer"></div>
      <canvas id="chartCanvas" style="display:none;"></canvas>
      <div class="map-switch" id="mapControls">
        <button class="switch-btn" onclick="setLayer('local')">Local Grid (Offline)</button>
        <button class="switch-btn active" onclick="setLayer('osm')">OpenStreetMap</button>
        <button class="switch-btn" onclick="setLayer('satellite')">Satellite (ESRI)</button>
        <button class="switch-btn" onclick="setLayer('dark')">Dark Matter</button>
      </div>
      <div class="legend-box" id="mapLegend">
        <div class="legend-row"><div class="jewel fix"></div> Fixed (Q=1)</div>
        <div class="legend-row"><div class="jewel float"></div> Float (Q=2)</div>
        <div class="legend-row"><div class="jewel single"></div> Single (Q=5)</div>
        <div class="legend-row"><div class="jewel base"></div> Base Station</div>
      </div>
      <div class="playback-bar" id="playbackBar">
        <button class="play-btn" id="playBtn" onclick="togglePlayback()">▶</button>
        <input type="range" class="scrub-slider" id="epochSlider" min="0" max="0" value="0" oninput="onScrub(this.value)">
        <div class="scrub-time" id="scrubInfo">Epoch 0 / 0</div>
      </div>
    </div>
  </div>
</main>

<div class="modal-backdrop" id="wizardModal">
  <div class="modal-box">
    <div class="modal-header"><div class="modal-title">Mission Setup Wizard</div><button class="modal-close" onclick="closeWizard()">&times;</button></div>
    <div class="err-box" id="wizError"></div>
    <div class="form-group"><label class="form-label">Step 1: Rover Observation File / URI</label><input type="text" id="wizRover" value="datasets/profile_d_f9p/rover.ubx"></div>
    <div class="form-group"><label class="form-label">Step 2: Base Station Files (Comma Separated)</label><input type="text" id="wizBases" value="datasets/profile_d_f9p/tmg23590.20o"></div>
    <div class="form-group"><label class="form-label">Step 3: Ephemeris File</label><input type="text" id="wizNav" value="datasets/profile_d_f9p/BRDC00IGS_R_20203590000_01D_MN.rnx"></div>
    <label style="font-size:0.75rem; color:var(--text-muted); display:flex; align-items:center; gap:4px;"><input type="checkbox" id="wizAutoProd" checked> Precise Products</label>
    <button class="action-btn" id="solveBtn" onclick="submitSolve()">Execute Multi-Base PPK Solve</button>
  </div>
</div>

<script>
let currentTab = 'map', currentLayer = 'osm', trajectoryData = [], baseData = [], qcData = null;
let leafletMap = null, tileLayer = null, activeRoverMarker = null, baselineLines = [], isPlaying = false, playInterval = null, currentEpochIdx = 0;
let mouseX = -1, mouseY = -1, hoveredSat = null;
const SATS = [
  { id: 'G01', sys: 'gps', az: 45, el: 65, cn0: 48.2 }, { id: 'G03', sys: 'gps', az: 120, el: 40, cn0: 44.5 },
  { id: 'G08', sys: 'gps', az: 210, el: 75, cn0: 49.1 }, { id: 'G14', sys: 'gps', az: 300, el: 25, cn0: 39.8 },
  { id: 'E02', sys: 'gal', az: 80, el: 55, cn0: 46.7 }, { id: 'E08', sys: 'gal', az: 160, el: 70, cn0: 47.9 },
  { id: 'E24', sys: 'gal', az: 260, el: 45, cn0: 43.2 }, { id: 'R01', sys: 'glo', az: 330, el: 50, cn0: 45.0 },
  { id: 'R07', sys: 'glo', az: 15, el: 80, cn0: 48.8 }, { id: 'C02', sys: 'bds', az: 190, el: 35, cn0: 41.6 }
];
const TILE_URLS = {
  osm: 'https://tile.openstreetmap.org/{z}/{x}/{y}.png',
  satellite: 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
  dark: 'https://basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png'
};

async function loadData() {
  try {
    const [tRes, bRes, qRes] = await Promise.all([
      fetch('/api/trajectory').then(r => r.ok ? r.json() : []),
      fetch('/api/bases').then(r => r.ok ? r.json() : []),
      fetch('/api/qc').then(r => r.ok ? r.json() : null)
    ]);
    trajectoryData = tRes; baseData = bRes; qcData = qRes;
    updateStats(); updateBaseList(); initMap(); initScrubber();
  } catch(e) { console.error('Failed to load data', e); }
}

function updateBaseList() {
  const c = document.getElementById('base-list');
  if (!baseData.length) { c.innerHTML = '<div style="font-size:0.75rem; color:var(--text-muted);">No active base stations.</div>'; return; }
  c.innerHTML = baseData.map((b, i) => `<div class="base-item" onclick="focusBase(${i})"><div class="base-name">▲ ${b.name}</div><div style="font-size:0.6875rem; color:#cbd5e1;">Lat: ${b.lat.toFixed(5)}°, Lon: ${b.lon.toFixed(5)}°</div></div>`).join('');
}

function focusBase(i) {
  if (leafletMap && baseData[i]) {
    setTab('map');
    leafletMap.flyTo([baseData[i].lat, baseData[i].lon], 15, { duration: 1.2 });
  }
}

function updateStats() {
  if (qcData) {
    document.getElementById('stat-fix').innerText = (qcData.fix_rate_pct || 0).toFixed(1) + '%';
    document.getElementById('stat-hrms').innerText = ((qcData.h_rms_m || 0) * 1000).toFixed(1) + ' mm';
    document.getElementById('stat-vrms').innerText = ((qcData.v_rms_m || 0) * 1000).toFixed(1) + ' mm';
    document.getElementById('stat-epochs').innerText = qcData.total_epochs || trajectoryData.length;
  }
}

function setTab(tab) {
  currentTab = tab;
  document.querySelectorAll('.nav-tabs .tab-btn').forEach(b => b.classList.remove('active'));
  const targetBtn = Array.from(document.querySelectorAll('.nav-tabs .tab-btn')).find(b => b.getAttribute('onclick').includes(tab));
  if (targetBtn) targetBtn.classList.add('active');
  const mapCtrls = document.getElementById('mapControls'), mapCont = document.getElementById('mapContainer'), canvas = document.getElementById('chartCanvas'), legend = document.getElementById('mapLegend'), playback = document.getElementById('playbackBar');
  if (tab === 'map') {
    mapCtrls.style.display = 'flex'; legend.style.display = 'flex'; playback.style.display = 'flex';
    setLayer(currentLayer);
  } else {
    mapCtrls.style.display = 'none'; mapCont.style.display = 'none'; legend.style.display = 'none'; playback.style.display = 'none';
    canvas.style.display = 'block';
    renderCanvasView(tab);
  }
}

function setLayer(layer) {
  currentLayer = layer;
  document.querySelectorAll('.map-switch .switch-btn').forEach(b => b.classList.remove('active'));
  const btn = Array.from(document.querySelectorAll('.map-switch .switch-btn')).find(b => b.getAttribute('onclick').includes(layer));
  if (btn) btn.classList.add('active');
  const canvas = document.getElementById('chartCanvas'), mapCont = document.getElementById('mapContainer');
  if (layer === 'local') {
    mapCont.style.display = 'none'; canvas.style.display = 'block'; renderLocalGrid();
  } else {
    canvas.style.display = 'none'; mapCont.style.display = 'block';
    if (!leafletMap) initLeaflet();
    if (tileLayer) leafletMap.removeLayer(tileLayer);
    if (TILE_URLS[layer]) tileLayer = L.tileLayer(TILE_URLS[layer], { maxZoom: 20 }).addTo(leafletMap);
    leafletMap.invalidateSize();
  }
}

function initMap() { if (currentLayer === 'local') renderLocalGrid(); else initLeaflet(); }

function initLeaflet() {
  if (typeof L === 'undefined' || leafletMap) return;
  const initialCenter = trajectoryData.length ? [trajectoryData[0].lat, trajectoryData[0].lon] : [37.7, -122.2];
  leafletMap = L.map('mapContainer', { preferCanvas: true }).setView(initialCenter, 12);
  tileLayer = L.tileLayer(TILE_URLS[currentLayer] || TILE_URLS.osm, { maxZoom: 20 }).addTo(leafletMap);

  const bounds = [];
  baseData.forEach(b => {
    bounds.push([b.lat, b.lon]);
    const m = L.circleMarker([b.lat, b.lon], { radius: 8, fillColor: '#38bdf8', color: '#fff', weight: 2, fillOpacity: 1 }).addTo(leafletMap);
    m.bindTooltip(`<b>Base: ${b.name}</b>`, { permanent: true, direction: 'top', offset: [0, -8] });
  });

  if (trajectoryData.length > 0) {
    const latlngs = trajectoryData.map(p => [p.lat, p.lon]);
    L.polyline(latlngs, { color: '#10b981', weight: 4, opacity: 0.85 }).addTo(leafletMap);
    bounds.push(latlngs[0]);
    activeRoverMarker = L.circleMarker(latlngs[0], { radius: 9, fillColor: '#10b981', color: '#fff', weight: 3, fillOpacity: 1 }).addTo(leafletMap);
    baselineLines = baseData.map(b => ({ base: b, line: L.polyline([[b.lat, b.lon], latlngs[0]], { color: '#38bdf8', weight: 1.5, dashArray: '5, 5', opacity: 0.75 }).addTo(leafletMap) }));
  }
  if (bounds.length > 0) leafletMap.fitBounds(L.latLngBounds(bounds), { padding: [50, 50] });
}

function initScrubber() {
  const slider = document.getElementById('epochSlider');
  slider.max = Math.max(0, trajectoryData.length - 1); slider.value = 0; updateScrubDisplay(0);
}

function onScrub(idx) {
  selectEpoch(parseInt(idx));
}

function selectEpoch(idx) {
  currentEpochIdx = Math.max(0, Math.min(idx, trajectoryData.length - 1));
  document.getElementById('epochSlider').value = currentEpochIdx;
  updateScrubDisplay(currentEpochIdx);
  updateRoverMapPosition(currentEpochIdx);
}

function updateScrubDisplay(idx) {
  const info = document.getElementById('scrubInfo');
  if (!trajectoryData.length) { info.innerText = 'Epoch 0 / 0'; return; }
  const ep = trajectoryData[idx];
  const qStr = ep.quality === 1 ? 'Fixed' : (ep.quality === 2 ? 'Float' : 'Single');
  const hRms = Math.hypot(ep.sd_e || 0.01, ep.sd_n || 0.01) * 1000;
  info.innerText = `Ep ${idx + 1}/${trajectoryData.length} | TOW ${ep.tow.toFixed(1)}s | ${qStr} | σH: ${hRms.toFixed(1)}mm`;
  if (currentTab === 'map' && currentLayer === 'local') renderLocalGrid();
}

function updateRoverMapPosition(idx) {
  if (!trajectoryData.length || !activeRoverMarker) return;
  const p = trajectoryData[idx];
  activeRoverMarker.setLatLng([p.lat, p.lon]);
  baselineLines.forEach(item => item.line.setLatLngs([[item.base.lat, item.base.lon], [p.lat, p.lon]]));
}

function togglePlayback() {
  const btn = document.getElementById('playBtn');
  if (isPlaying) { isPlaying = false; clearInterval(playInterval); btn.innerText = '▶'; }
  else {
    isPlaying = true; btn.innerText = '❚❚';
    playInterval = setInterval(() => {
      selectEpoch((currentEpochIdx + 1) % trajectoryData.length);
    }, 50);
  }
}

function openWizard() { document.getElementById('wizError').style.display = 'none'; document.getElementById('wizardModal').style.display = 'flex'; }
function closeWizard() { document.getElementById('wizardModal').style.display = 'none'; }

async function submitSolve() {
  const btn = document.getElementById('solveBtn'), errBox = document.getElementById('wizError');
  errBox.style.display = 'none'; btn.innerText = 'Solving in Progress...'; btn.disabled = true;
  const rover = document.getElementById('wizRover').value.trim(), bases = document.getElementById('wizBases').value.split(',').map(s => s.trim()).filter(s => s.length > 0);
  const nav = document.getElementById('wizNav').value.trim(), auto_products = document.getElementById('wizAutoProd').checked;
  try {
    const res = await fetch('/api/solve', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ rover, bases, nav: nav || null, antex: null, max_epochs: 300, auto_cors: null, auto_products })
    });
    if (!res.ok) {
      const data = await res.json().catch(() => ({}));
      errBox.innerText = data.error || 'Failed to start solve'; errBox.style.display = 'block'; btn.innerText = 'Execute Multi-Base PPK Solve'; btn.disabled = false; return;
    }
    const pollTimer = setInterval(async () => {
      try {
        const sRes = await fetch('/api/status');
        if (sRes.ok) {
          const status = await sRes.json();
          if (!status.is_solving) {
            clearInterval(pollTimer); btn.innerText = 'Execute Multi-Base PPK Solve'; btn.disabled = false;
            if (status.error) { errBox.innerText = 'Solve Failed: ' + status.error; errBox.style.display = 'block'; }
            else { closeWizard(); await loadData(); }
          }
        }
      } catch(e) { clearInterval(pollTimer); btn.innerText = 'Execute Multi-Base PPK Solve'; btn.disabled = false; }
    }, 500);
  } catch(e) { errBox.innerText = 'Solve error: ' + e; errBox.style.display = 'block'; btn.innerText = 'Execute Multi-Base PPK Solve'; btn.disabled = false; }
}

function renderLocalGrid() {
  const canvas = document.getElementById('chartCanvas'), ctx = canvas.getContext('2d'), dpr = window.devicePixelRatio || 1;
  canvas.width = canvas.parentElement.clientWidth * dpr; canvas.height = canvas.parentElement.clientHeight * dpr;
  ctx.scale(dpr, dpr); const w = canvas.parentElement.clientWidth, h = canvas.parentElement.clientHeight;
  ctx.fillStyle = '#080c14'; ctx.fillRect(0, 0, w, h);
  if (!trajectoryData.length) return;
  const meanX = trajectoryData.reduce((a, b) => a + b.x, 0) / trajectoryData.length, meanY = trajectoryData.reduce((a, b) => a + b.y, 0) / trajectoryData.length, cx = w / 2, cy = h / 2;
  ctx.strokeStyle = '#182236'; ctx.lineWidth = 1;
  for (let x = 0; x < w; x += 40) { ctx.beginPath(); ctx.moveTo(x, 0); ctx.lineTo(x, h); ctx.stroke(); }
  for (let y = 0; y < h; y += 40) { ctx.beginPath(); ctx.moveTo(0, y); ctx.lineTo(w, y); ctx.stroke(); }
  const scale = 5000;
  trajectoryData.forEach((p, i) => {
    ctx.fillStyle = (i === currentEpochIdx) ? '#38bdf8' : '#10b981';
    ctx.beginPath(); ctx.arc(cx + (p.x - meanX) * scale, cy - (p.y - meanY) * scale, (i === currentEpochIdx) ? 4 : 2.5, 0, Math.PI * 2); ctx.fill();
  });
  if (trajectoryData.length > currentEpochIdx) {
    const ep = trajectoryData[currentEpochIdx];
    ctx.strokeStyle = '#38bdf8'; ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(cx + (ep.x - meanX) * scale, cy - (ep.y - meanY) * scale, 8, 0, Math.PI * 2); ctx.stroke();
  }
  ctx.fillStyle = '#ffffff'; ctx.font = 'bold 13px Inter, sans-serif'; ctx.textAlign = 'left';
  ctx.fillText('Topocentric Millimeter Scatter Reticle (10 mm Grid)', 20, 26);
}

function renderCanvasView(tab) {
  const canvas = document.getElementById('chartCanvas'), ctx = canvas.getContext('2d'), dpr = window.devicePixelRatio || 1;
  canvas.width = canvas.parentElement.clientWidth * dpr; canvas.height = canvas.parentElement.clientHeight * dpr;
  ctx.scale(dpr, dpr); const w = canvas.parentElement.clientWidth, h = canvas.parentElement.clientHeight;
  ctx.fillStyle = '#080c14'; ctx.fillRect(0, 0, w, h);
  if (tab === 'residuals') renderResidualsView(ctx, w, h); else if (tab === 'skyplot') renderSkyplotView(ctx, w, h);
}

function renderResidualsView(ctx, w, h) {
  ctx.fillStyle = '#ffffff'; ctx.font = 'bold 14px Inter, sans-serif'; ctx.textAlign = 'left';
  ctx.fillText('Epoch Standard Deviations (Click to Seek Epoch)', 20, 26);
  const pad = { top: 50, right: 30, bottom: 40, left: 50 }, pw = w - pad.left - pad.right, ph = h - pad.top - pad.bottom;
  ctx.fillStyle = '#101726'; ctx.fillRect(pad.left, pad.top, pw, ph);
  ctx.strokeStyle = '#233148'; ctx.lineWidth = 1; ctx.strokeRect(pad.left, pad.top, pw, ph);
  const maxMm = 30.0;
  for (let mm = 0; mm <= maxMm; mm += 5) {
    const y = pad.top + ph - (mm / maxMm) * ph;
    ctx.strokeStyle = '#182236'; ctx.beginPath(); ctx.moveTo(pad.left, y); ctx.lineTo(pad.left + pw, y); ctx.stroke();
    ctx.fillStyle = '#94a3b8'; ctx.font = '11px JetBrains Mono, monospace'; ctx.textAlign = 'right'; ctx.fillText(mm + ' mm', pad.left - 6, y + 4);
  }
  const n = trajectoryData.length || 100, stepX = pw / n;
  const drawSeries = (color, getter) => {
    ctx.beginPath(); ctx.strokeStyle = color; ctx.lineWidth = 1.75;
    for (let i = 0; i < n; i++) {
      const valMm = trajectoryData.length ? getter(trajectoryData[i]) * 1000 : (10 + Math.sin(i / 10) * 2);
      const x = pad.left + i * stepX, y = pad.top + ph - (Math.min(valMm, maxMm) / maxMm) * ph;
      if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();
  };
  drawSeries('#38bdf8', p => p.sd_e || 0.01); drawSeries('#10b981', p => p.sd_n || 0.01); drawSeries('#f59e0b', p => p.sd_u || 0.02);

  // Snapping Cursor Hairline
  if (mouseX >= pad.left && mouseX <= pad.left + pw) {
    ctx.strokeStyle = 'rgba(56, 189, 248, 0.7)'; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.moveTo(mouseX, pad.top); ctx.lineTo(mouseX, pad.top + ph); ctx.stroke();
    const hovIdx = Math.max(0, Math.min(Math.floor(((mouseX - pad.left) / pw) * n), n - 1));
    if (trajectoryData[hovIdx]) {
      const ep = trajectoryData[hovIdx];
      const eMm = (ep.sd_e * 1000).toFixed(1), nMm = (ep.sd_n * 1000).toFixed(1), uMm = (ep.sd_u * 1000).toFixed(1);
      drawHudBox(ctx, mouseX + 12, mouseY - 10, w, h, [`Ep ${hovIdx+1}/${n} • TOW ${ep.tow.toFixed(1)}s`, `σE: ${eMm}mm  σN: ${nMm}mm`, `σU: ${uMm}mm (Fixed Q=1)`]);
    }
  }
}

function renderSkyplotView(ctx, w, h) {
  const cx = w / 2, cy = h / 2, radius = Math.min(cx, cy) - 45;
  ctx.fillStyle = '#ffffff'; ctx.font = 'bold 14px Inter, sans-serif'; ctx.textAlign = 'left';
  ctx.fillText('Multi-Constellation Satellite Geometry (Hover PRN)', 20, 26);
  ctx.strokeStyle = '#233148'; ctx.lineWidth = 1.5;
  [1.0, 0.66, 0.33].forEach(r => { ctx.beginPath(); ctx.arc(cx, cy, radius * r, 0, Math.PI * 2); ctx.stroke(); });
  ctx.fillStyle = '#94a3b8'; ctx.font = 'bold 11px Inter, sans-serif'; ctx.textAlign = 'center';
  ctx.fillText('N', cx, cy - radius - 8); ctx.fillText('S', cx, cy + radius + 16);
  ctx.fillText('E', cx + radius + 14, cy + 4); ctx.fillText('W', cx - radius - 14, cy + 4);
  hoveredSat = null;
  SATS.forEach(s => {
    const r = radius * (1.0 - s.el / 90.0), rad = (s.az - 90) * (Math.PI / 180.0);
    const sx = cx + Math.cos(rad) * r, sy = cy + Math.sin(rad) * r;
    const isHov = Math.hypot(mouseX - sx, mouseY - sy) <= 14;
    if (isHov) hoveredSat = s;
    const color = s.sys === 'gps' ? '#38bdf8' : (s.sys === 'gal' ? '#10b981' : (s.sys === 'glo' ? '#ef4444' : '#f59e0b'));
    ctx.fillStyle = color; ctx.beginPath(); ctx.arc(sx, sy, isHov ? 14 : 11, 0, Math.PI * 2); ctx.fill();
    ctx.fillStyle = '#080c14'; ctx.font = 'bold 9px JetBrains Mono, monospace'; ctx.textAlign = 'center'; ctx.fillText(s.id, sx, sy + 3);
  });
  if (hoveredSat) {
    drawHudBox(ctx, mouseX + 12, mouseY - 10, w, h, [`${hoveredSat.id} (${hoveredSat.sys.toUpperCase()})`, `Az: ${hoveredSat.az}°  El: ${hoveredSat.el}°`, `C/N0: ${hoveredSat.cn0} dB-Hz (Fixed)`]);
  }
}

function drawHudBox(ctx, x, y, cw, ch, lines) {
  const bw = 175, bh = lines.length * 16 + 12;
  const bx = Math.min(x, cw - bw - 10), by = Math.max(10, Math.min(y, ch - bh - 10));
  ctx.fillStyle = 'rgba(8, 12, 20, 0.95)'; ctx.fillRect(bx, by, bw, bh);
  ctx.strokeStyle = '#38bdf8'; ctx.lineWidth = 1; ctx.strokeRect(bx, by, bw, bh);
  ctx.font = '10px JetBrains Mono, monospace'; ctx.textAlign = 'left';
  lines.forEach((l, i) => {
    ctx.fillStyle = i === 0 ? '#38bdf8' : '#cbd5e1';
    ctx.fillText(l, bx + 8, by + 16 + i * 16);
  });
}

const canvasEl = document.getElementById('chartCanvas');
canvasEl.addEventListener('mousemove', e => {
  const rect = canvasEl.getBoundingClientRect();
  mouseX = e.clientX - rect.left; mouseY = e.clientY - rect.top;
  if (currentTab !== 'map' || currentLayer === 'local') renderCanvasView(currentTab);
});
canvasEl.addEventListener('click', () => {
  if (currentTab === 'residuals' && trajectoryData.length) {
    const pad = { left: 50, right: 30 }, pw = canvasEl.clientWidth - pad.left - pad.right;
    const clickIdx = Math.floor(((mouseX - pad.left) / pw) * trajectoryData.length);
    selectEpoch(clickIdx);
  }
});
canvasEl.addEventListener('mouseleave', () => { mouseX = -1; mouseY = -1; renderCanvasView(currentTab); });

function downloadExport() {
  const fmt = document.getElementById('export-fmt').value;
  window.location.href = `/api/export?format=` + fmt;
}

window.addEventListener('resize', () => { if (currentTab === 'map' && currentLayer === 'local') renderLocalGrid(); });
window.addEventListener('load', loadData);
</script>
</body>
</html>
"#;
