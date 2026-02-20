"""
AD-HTC Fuel-Enhanced Gas Cycle — Process Flow Analysis Dashboard
================================================================
A Dash-based interactive web application for analyzing an Anaerobic Digestion
and Hydrothermal Carbonization (AD-HTC) fuel-enhanced gas power cycle.

Features:
  - Animated SVG schematic of the full AD-HTC process
  - Dual biomass tank input controls (Tank A: moisture-lean, Tank B: moisture-rich)
  - Brayton (gas turbine) cycle thermodynamic analysis
  - HTC steam (Rankine-like) cycle thermodynamic analysis
  - AD biogas yield estimation
  - h-s diagram for the HTC steam cycle
  - T-Ḣ (temperature vs enthalpy rate) chart for the gas cycle
  - Comprehensive report with state-point tables and energy balance

Run:
    pip install -r requirements.txt
    python app.py

Then open http://127.0.0.1:8050 in your browser.
"""

import math
import io
import base64

import dash
from dash import html, dcc, Input, Output, State, dash_table, callback

import matplotlib
matplotlib.use("Agg")  # non-interactive backend — must be set before importing pyplot
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# ═══════════════════════════════════════════════════════════════════════════════
# THERMODYNAMIC CONSTANTS
# ═══════════════════════════════════════════════════════════════════════════════

CP_AIR = 1.005       # kJ/(kg·K) — specific heat of air at constant pressure
GAMMA_AIR = 1.4      # ratio of specific heats for air
R_AIR = 0.287        # kJ/(kg·K) — gas constant for air
CP_STEAM = 2.01      # kJ/(kg·K) — approx for superheated steam
HFG_WATER = 2257.0   # kJ/kg — latent heat of vaporization
CP_WATER = 4.186     # kJ/(kg·K) — specific heat of water
BIOGAS_LHV = 22.0    # MJ/m³ — lower heating value of biogas
HYDROCHAR_HHV = 25.0 # MJ/kg — higher heating value of hydrochar


# ═══════════════════════════════════════════════════════════════════════════════
# THERMODYNAMIC CALCULATIONS
# ═══════════════════════════════════════════════════════════════════════════════

def brayton_cycle(T1_C, rp, T3_C, eta_c, eta_t):
    """
    Calculate Brayton (gas turbine) cycle state points and performance.

    Parameters:
        T1_C  : Ambient / compressor inlet temperature (°C)
        rp    : Compressor pressure ratio
        T3_C  : Turbine inlet temperature (°C)
        eta_c : Isentropic efficiency of compressor (fraction, e.g. 0.85)
        eta_t : Isentropic efficiency of turbine (fraction, e.g. 0.90)

    Returns:
        states  : dict with T, h, s arrays and labels for each state point
        metrics : dict with work, heat, and efficiency values
    """
    T1 = T1_C + 273.15
    T3 = T3_C + 273.15

    # Isentropic compression
    T2s = T1 * rp ** ((GAMMA_AIR - 1) / GAMMA_AIR)
    T2 = T1 + (T2s - T1) / eta_c

    # Isentropic expansion
    T4s = T3 / rp ** ((GAMMA_AIR - 1) / GAMMA_AIR)
    T4 = T3 - (T3 - T4s) * eta_t

    # Work & heat (per kg of air)
    w_comp = CP_AIR * (T2 - T1)
    w_turb = CP_AIR * (T3 - T4)
    w_net = w_turb - w_comp
    q_in = CP_AIR * (T3 - T2)
    q_out = CP_AIR * (T4 - T1)
    eta_th = (w_net / q_in * 100) if q_in > 0 else 0
    bwr = (w_comp / w_turb * 100) if w_turb > 0 else 0

    # Entropy (relative to state 1 = 0)
    s1 = 0.0
    s2 = CP_AIR * math.log(T2 / T1) - R_AIR * math.log(rp)
    s3 = s2 + CP_AIR * math.log(T3 / T2)
    s4 = s3 + CP_AIR * math.log(T4 / T3) + R_AIR * math.log(rp)

    # Enthalpy (relative to state 1 = 0)
    h1 = 0.0
    h2 = CP_AIR * (T2 - T1)
    h3 = CP_AIR * (T3 - T1)
    h4 = CP_AIR * (T4 - T1)

    states = dict(
        T=[T1 - 273.15, T2 - 273.15, T3 - 273.15, T4 - 273.15],
        h=[h1, h2, h3, h4],
        s=[s1, s2, s3, s4],
        labels=["1 – Comp Inlet", "2 – Comp Outlet",
                "3 – Turb Inlet", "4 – Turb Outlet"],
    )

    metrics = dict(
        w_comp=round(w_comp, 2),
        w_turb=round(w_turb, 2),
        w_net=round(w_net, 2),
        q_in=round(q_in, 2),
        q_out=round(q_out, 2),
        eta_th=round(eta_th, 2),
        bwr=round(bwr, 2),
    )

    return states, metrics


def htc_steam_cycle(T_reactor=200, P_reactor=20):
    """
    Simplified HTC steam (Rankine-like) cycle analysis.

    Parameters:
        T_reactor : HTC reactor temperature (°C)
        P_reactor : HTC reactor pressure (bar)

    Returns:
        states  : dict with T, h, s arrays for each state point
        metrics : dict with work, heat, and efficiency values
    """
    T1 = 45.0          # condenser temperature (°C)
    T3 = T_reactor + 50  # superheat above reactor temp
    P_low = 0.1        # condenser pressure (bar)
    P_high = P_reactor

    # Approximate enthalpies
    h1 = CP_WATER * T1
    h2 = h1 + 0.001 * (P_high - P_low) * 100   # pump work (small)
    h3 = CP_WATER * 100 + HFG_WATER + CP_STEAM * (T3 - 100)

    # Turbine with 85 % isentropic efficiency
    eta_st = 0.85
    h4s = h1 + 0.88 * HFG_WATER
    h4 = h3 - eta_st * (h3 - h4s)

    # Approximate entropies
    s1 = CP_WATER * math.log((T1 + 273.15) / 273.15)
    s2 = s1  # isentropic pump
    s3 = s1 + (h3 - h2) / (T3 + 273.15)
    s4 = s3 + 0.15  # entropy increase due to irreversibilities

    states = dict(
        T=[T1, T1 + 2, T3, T1 + 20],
        h=[round(h1, 1), round(h2, 1), round(h3, 1), round(h4, 1)],
        s=[round(s1, 4), round(s2, 4), round(s3, 4), round(s4, 4)],
        labels=["1 – Pump Inlet", "2 – Pump Outlet",
                "3 – Boiler Outlet", "4 – Turb Outlet"],
    )

    w_pump = h2 - h1
    w_turb_st = h3 - h4
    q_in_st = h3 - h2
    w_net_st = w_turb_st - w_pump
    eta = (w_net_st / q_in_st * 100) if q_in_st > 0 else 0

    metrics = dict(
        w_pump=round(w_pump, 2),
        w_turb=round(w_turb_st, 2),
        w_net=round(w_net_st, 2),
        q_in=round(q_in_st, 2),
        eta=round(eta, 2),
    )

    return states, metrics


def ad_biogas_yield(feed_rate, moisture_pct, vs_fraction=0.80):
    """
    Estimate biogas yield from anaerobic digestion of moisture-rich biomass.
    """
    dry_mass = feed_rate * (1 - moisture_pct / 100)
    volatile_solids = dry_mass * vs_fraction
    biogas_m3_h = volatile_solids * 0.40       # m³/h typical yield
    methane_m3_h = biogas_m3_h * 0.60          # 60 % methane
    biogas_energy = biogas_m3_h * BIOGAS_LHV   # MJ/h

    return dict(
        dry_mass=round(dry_mass, 2),
        volatile_solids=round(volatile_solids, 2),
        biogas_m3_h=round(biogas_m3_h, 2),
        methane_m3_h=round(methane_m3_h, 2),
        biogas_energy_MJ_h=round(biogas_energy, 2),
    )


def htc_process(feed_rate, moisture_pct, T_reactor=200):
    """
    Estimate HTC (Hydrothermal Carbonization) products from moisture-lean biomass.
    """
    dry_mass = feed_rate * (1 - moisture_pct / 100)
    hydrochar_yield = dry_mass * 0.60
    process_water = feed_rate - hydrochar_yield
    hydrochar_energy = hydrochar_yield * HYDROCHAR_HHV   # MJ/h
    energy_required = feed_rate * CP_WATER * (T_reactor - 25) / 1000  # MJ/h

    return dict(
        dry_mass=round(dry_mass, 2),
        hydrochar_kg_h=round(hydrochar_yield, 2),
        process_water_kg_h=round(process_water, 2),
        hydrochar_energy_MJ_h=round(hydrochar_energy, 2),
        energy_required_MJ_h=round(energy_required, 2),
    )


# ═══════════════════════════════════════════════════════════════════════════════
# ANIMATED SVG SCHEMATIC (inline HTML)
# ═══════════════════════════════════════════════════════════════════════════════

def build_schematic_html():
    """Return a self-contained HTML document with the animated AD-HTC SVG."""
    return """<!DOCTYPE html>
<html>
<head>
<style>
  * { margin:0; padding:0; box-sizing:border-box; }
  body {
    background: #0d0d20;
    display: flex;
    align-items: center;
    justify-content: center;
    min-height: 100vh;
    overflow: hidden;
    font-family: 'Inter', 'Segoe UI', sans-serif;
  }
  svg { width: 100%; height: auto; max-height: 95vh; }

  /* ── Animated flow dashes ── */
  .flow-pipe {
    fill: none;
    stroke-width: 3;
    stroke-dasharray: 12 8;
    animation: pipe-flow 1.2s linear infinite;
  }
  @keyframes pipe-flow {
    to { stroke-dashoffset: -40; }
  }
  .flow-biomass { stroke: #8B6914; }
  .flow-water   { stroke: #00bfff; }
  .flow-biogas  { stroke: #44dd88; }
  .flow-steam   { stroke: #ff6644; }
  .flow-air     { stroke: #aabbff; }
  .flow-exhaust { stroke: #ff8844; opacity:0.7; }

  /* ── Component boxes ── */
  .comp-box {
    fill: rgba(20,20,50,0.85);
    stroke-width: 2;
    rx: 8; ry: 8;
    transition: all 0.3s;
  }
  .comp-box:hover { filter: brightness(1.3); }

  .comp-htc    { stroke: #00bfff; }
  .comp-ad     { stroke: #44dd88; }
  .comp-turb   { stroke: #ff8833; }
  .comp-combust{ stroke: #ff4466; }
  .comp-boiler { stroke: #ff6644; }
  .comp-homog  { stroke: #8B6914; }
  .comp-biogas { stroke: #44dd88; }
  .comp-feed   { stroke: #b8860b; }

  /* Tank styles */
  .tank {
    fill: rgba(20,20,50,0.9);
    stroke-width: 2.5;
    rx: 12; ry: 12;
  }
  .tank-a { stroke: #00d4ff; }
  .tank-b { stroke: #00e88f; }

  /* ── Labels ── */
  .label {
    fill: #e0e0f0;
    font-size: 11px;
    font-weight: 600;
    text-anchor: middle;
    dominant-baseline: central;
    pointer-events: none;
  }
  .label-sm {
    fill: #9999bb;
    font-size: 9px;
    font-weight: 400;
    text-anchor: middle;
    dominant-baseline: central;
  }
  .label-title {
    fill: #ffffff;
    font-size: 16px;
    font-weight: 800;
    text-anchor: middle;
    letter-spacing: 1px;
  }
  .label-tank {
    fill: #ffffff;
    font-size: 13px;
    font-weight: 700;
    text-anchor: middle;
    dominant-baseline: central;
  }
  .io-label {
    fill: #8888aa;
    font-size: 10px;
    font-weight: 500;
    text-anchor: middle;
  }

  /* ── Pulsing animation ── */
  .pulse {
    animation: comp-pulse 2.5s ease-in-out infinite;
  }
  @keyframes comp-pulse {
    0%, 100% { opacity: 0.85; }
    50%      { opacity: 1; filter: drop-shadow(0 0 8px rgba(0,212,255,0.3)); }
  }

  /* ── Rotating turbine ── */
  .rotate-cw {
    transform-origin: center;
    animation: spin-cw 3s linear infinite;
  }
  @keyframes spin-cw {
    to { transform: rotate(360deg); }
  }

  /* Arrow markers */
  .arrow-biomass { fill: #8B6914; }
  .arrow-water   { fill: #00bfff; }
  .arrow-biogas  { fill: #44dd88; }
  .arrow-steam   { fill: #ff6644; }
  .arrow-air     { fill: #aabbff; }

  /* Tank fill animation */
  .tank-fill-a {
    fill: rgba(0,212,255,0.15);
    animation: fill-pulse 3s ease-in-out infinite;
  }
  .tank-fill-b {
    fill: rgba(0,232,143,0.15);
    animation: fill-pulse 3s ease-in-out infinite 1s;
  }
  @keyframes fill-pulse {
    0%, 100% { opacity: 0.4; }
    50% { opacity: 0.8; }
  }
</style>
</head>
<body>
<svg viewBox="0 0 1200 640" xmlns="http://www.w3.org/2000/svg">
  <!-- ── Defs: arrow markers ── -->
  <defs>
    <marker id="arr-bio" markerWidth="8" markerHeight="6" refX="8" refY="3" orient="auto">
      <path d="M0,0 L8,3 L0,6 Z" class="arrow-biomass"/>
    </marker>
    <marker id="arr-water" markerWidth="8" markerHeight="6" refX="8" refY="3" orient="auto">
      <path d="M0,0 L8,3 L0,6 Z" class="arrow-water"/>
    </marker>
    <marker id="arr-gas" markerWidth="8" markerHeight="6" refX="8" refY="3" orient="auto">
      <path d="M0,0 L8,3 L0,6 Z" class="arrow-biogas"/>
    </marker>
    <marker id="arr-steam" markerWidth="8" markerHeight="6" refX="8" refY="3" orient="auto">
      <path d="M0,0 L8,3 L0,6 Z" class="arrow-steam"/>
    </marker>
    <marker id="arr-air" markerWidth="8" markerHeight="6" refX="8" refY="3" orient="auto">
      <path d="M0,0 L8,3 L0,6 Z" class="arrow-air"/>
    </marker>

    <!-- glow filter -->
    <filter id="glow">
      <feGaussianBlur stdDeviation="3" result="blur"/>
      <feMerge><feMergeNode in="blur"/><feMergeNode in="SourceGraphic"/></feMerge>
    </filter>
  </defs>

  <!-- ── Title ── -->
  <text x="600" y="32" class="label-title">AD-HTC Fuel-Enhanced Gas Power Cycle Schematic</text>

  <!-- ════════════════════════════════════════════════════ -->
  <!-- BIOMASS FEEDSTOCK (input) -->
  <!-- ════════════════════════════════════════════════════ -->
  <text x="62" y="135" class="io-label">Biomass</text>
  <text x="62" y="148" class="io-label">Feedstock</text>
  <line x1="100" y1="140" x2="148" y2="140" class="flow-pipe flow-biomass" marker-end="url(#arr-bio)"/>

  <!-- HOMOGENIZER -->
  <rect x="150" y="110" width="140" height="60" class="comp-box comp-homog pulse"/>
  <text x="220" y="136" class="label">Biomass Feedstock</text>
  <text x="220" y="150" class="label">Homogenizer</text>

  <!-- Split from homogenizer -->
  <!-- Down to Tank A (moisture-lean) -->
  <line x1="220" y1="170" x2="220" y2="208" class="flow-pipe flow-biomass" marker-end="url(#arr-bio)"/>
  <!-- Right to Tank B (moisture-rich) -->
  <line x1="290" y1="140" x2="535" y2="140" class="flow-pipe flow-biomass"/>
  <line x1="535" y1="140" x2="535" y2="148" class="flow-pipe flow-biomass" marker-end="url(#arr-bio)"/>

  <!-- ════════════════════════════════════════════════════ -->
  <!-- TANK A — Moisture-lean Biomass -->
  <!-- ════════════════════════════════════════════════════ -->
  <rect x="155" y="210" width="130" height="80" class="tank tank-a"/>
  <rect x="160" y="245" width="120" height="40" class="tank-fill-a" rx="6" ry="6"/>
  <text x="220" y="237" class="label-tank" fill="#00d4ff">TANK A</text>
  <text x="220" y="275" class="label-sm">Moisture-lean</text>
  <text x="220" y="286" class="label-sm">Biomass</text>

  <!-- Tank A → HTC Reactor -->
  <line x1="220" y1="290" x2="220" y2="338" class="flow-pipe flow-biomass" marker-end="url(#arr-bio)"/>

  <!-- ════════════════════════════════════════════════════ -->
  <!-- TANK B — Moisture-rich Biomass -->
  <!-- ════════════════════════════════════════════════════ -->
  <rect x="470" y="150" width="130" height="80" class="tank tank-b"/>
  <rect x="475" y="185" width="120" height="40" class="tank-fill-b" rx="6" ry="6"/>
  <text x="535" y="177" class="label-tank" fill="#00e88f">TANK B</text>
  <text x="535" y="213" class="label-sm">Moisture-rich</text>
  <text x="535" y="224" class="label-sm">Biomass</text>

  <!-- Tank B → AD -->
  <line x1="600" y1="190" x2="668" y2="190" class="flow-pipe flow-biogas" marker-end="url(#arr-gas)"/>

  <!-- ════════════════════════════════════════════════════ -->
  <!-- HTC REACTOR -->
  <!-- ════════════════════════════════════════════════════ -->
  <rect x="155" y="340" width="130" height="55" class="comp-box comp-htc pulse"/>
  <text x="220" y="364" class="label">HTC Reactor</text>
  <text x="220" y="378" class="label-sm">~200°C, 20 bar</text>

  <!-- HTC Reactor ↔ Boiler (steam cycle) -->
  <line x1="220" y1="340" x2="220" y2="320" class="flow-pipe flow-steam"/>
  <line x1="220" y1="320" x2="370" y2="320" class="flow-pipe flow-steam"/>
  <line x1="370" y1="320" x2="370" y2="338" class="flow-pipe flow-steam" marker-end="url(#arr-steam)"/>

  <!-- Boiler return to HTC -->
  <line x1="370" y1="395" x2="370" y2="410" class="flow-pipe flow-water"/>
  <line x1="370" y1="410" x2="220" y2="410" class="flow-pipe flow-water"/>
  <line x1="220" y1="410" x2="220" y2="398" class="flow-pipe flow-water" marker-end="url(#arr-water)"/>

  <!-- HTC Steam Cycle label -->
  <text x="295" y="308" class="label-sm" fill="#ff6644">HTC Steam Cycle</text>

  <!-- ════════════════════════════════════════════════════ -->
  <!-- BOILER -->
  <!-- ════════════════════════════════════════════════════ -->
  <rect x="320" y="340" width="100" height="55" class="comp-box comp-boiler pulse"/>
  <text x="370" y="364" class="label">Boiler</text>

  <!-- ════════════════════════════════════════════════════ -->
  <!-- Volatile Matters & Feedstock Waste (output from HTC) -->
  <!-- ════════════════════════════════════════════════════ -->
  <line x1="155" y1="380" x2="100" y2="380" class="flow-pipe flow-biomass"/>
  <line x1="100" y1="380" x2="100" y2="430" class="flow-pipe flow-biomass" marker-end="url(#arr-bio)"/>
  <text x="100" y="445" class="io-label">Volatile Matters</text>
  <text x="100" y="458" class="io-label">& Feedstock Waste</text>

  <!-- Process water from HTC to Enhanced Biogas -->
  <line x1="285" y1="367" x2="310" y2="367" class="flow-pipe flow-water"/>
  <line x1="310" y1="367" x2="310" y2="250" class="flow-pipe flow-water"/>
  <line x1="310" y1="250" x2="418" y2="250" class="flow-pipe flow-water" marker-end="url(#arr-water)"/>

  <!-- ════════════════════════════════════════════════════ -->
  <!-- ANAEROBIC DIGESTER (AD) -->
  <!-- ════════════════════════════════════════════════════ -->
  <rect x="670" y="155" width="120" height="70" class="comp-box comp-ad pulse"/>
  <text x="730" y="183" class="label">Anaerobic</text>
  <text x="730" y="198" class="label">Digester (AD)</text>

  <!-- AD → Enhanced Biogas Collection -->
  <line x1="730" y1="225" x2="730" y2="248" class="flow-pipe flow-biogas"/>
  <line x1="730" y1="248" x2="562" y2="248" class="flow-pipe flow-biogas" marker-end="url(#arr-gas)"/>

  <!-- ════════════════════════════════════════════════════ -->
  <!-- ENHANCED BIOGAS COLLECTION -->
  <!-- ════════════════════════════════════════════════════ -->
  <rect x="420" y="232" width="140" height="45" class="comp-box comp-biogas pulse"/>
  <text x="490" y="250" class="label">Enhanced Biogas</text>
  <text x="490" y="263" class="label">Collection</text>

  <!-- Biogas → Combustion Chamber  -->
  <line x1="490" y1="277" x2="490" y2="430" class="flow-pipe flow-biogas"/>
  <line x1="490" y1="430" x2="540" y2="430" class="flow-pipe flow-biogas"/>
  <line x1="540" y1="430" x2="540" y2="478" class="flow-pipe flow-biogas" marker-end="url(#arr-gas)"/>

  <!-- Biogas → Distribution -->
  <line x1="560" y1="254" x2="730" y2="254" class="flow-pipe flow-biogas"/>
  <line x1="730" y1="254" x2="730" y2="268" class="flow-pipe flow-biogas" marker-end="url(#arr-gas)"/>
  <text x="730" y="285" class="io-label">Biogas Distribution</text>
  <text x="730" y="298" class="io-label">to Building Envelopes</text>

  <!-- ════════════════════════════════════════════════════ -->
  <!-- GAS TURBINE CYCLE -->
  <!-- ════════════════════════════════════════════════════ -->

  <!-- COMPRESSOR -->
  <rect x="320" y="480" width="120" height="55" class="comp-box comp-turb pulse"/>
  <text x="380" y="504" class="label">Compressor</text>

  <!-- Air → Compressor -->
  <text x="380" y="575" class="io-label">Air Intake</text>
  <line x1="380" y1="562" x2="380" y2="538" class="flow-pipe flow-air" marker-end="url(#arr-air)"/>

  <!-- Compressor → Combustion Chamber -->
  <line x1="440" y1="507" x2="478" y2="507" class="flow-pipe flow-air" marker-end="url(#arr-air)"/>

  <!-- COMBUSTION CHAMBER -->
  <rect x="480" y="480" width="130" height="55" class="comp-box comp-combust pulse"/>
  <text x="545" y="500" class="label">Biogas Combustion</text>
  <text x="545" y="515" class="label">Chamber</text>

  <!-- Combustion Chamber → Turbine -->
  <line x1="610" y1="507" x2="648" y2="507" class="flow-pipe flow-steam" marker-end="url(#arr-steam)"/>

  <!-- TURBINE -->
  <rect x="650" y="480" width="120" height="55" class="comp-box comp-turb pulse"/>
  <text x="710" y="504" class="label">Gas Turbine</text>
  <!-- Rotating element inside turbine -->
  <g class="rotate-cw" style="transform-origin: 710px 515px;">
    <line x1="695" y1="515" x2="725" y2="515" stroke="#ff8833" stroke-width="2"/>
    <line x1="710" y1="508" x2="710" y2="522" stroke="#ff8833" stroke-width="2"/>
  </g>

  <!-- Turbine → Exhaust -->
  <line x1="770" y1="507" x2="820" y2="507" class="flow-pipe flow-exhaust" marker-end="url(#arr-steam)"/>
  <line x1="820" y1="507" x2="820" y2="545" class="flow-pipe flow-exhaust"/>
  <text x="820" y="565" class="io-label">Exhaust</text>
  <text x="820" y="578" class="io-label">Gases</text>

  <!-- Turbine shaft → Compressor (work) -->
  <line x1="380" y1="540" x2="380" y2="555" stroke="#555577" stroke-width="1" stroke-dasharray="3 3"/>
  <line x1="380" y1="555" x2="710" y2="555" stroke="#555577" stroke-width="1" stroke-dasharray="3 3"/>
  <line x1="710" y1="540" x2="710" y2="555" stroke="#555577" stroke-width="1" stroke-dasharray="3 3"/>
  <text x="545" y="567" class="label-sm" fill="#555577">⟵ Shaft Work ⟶</text>

  <!-- ════════════════════════════════════════════════════ -->
  <!-- Power output indicator -->
  <!-- ════════════════════════════════════════════════════ -->
  <rect x="850" y="470" width="120" height="50" class="comp-box" stroke="#8855ff" filter="url(#glow)"/>
  <text x="910" y="490" class="label" fill="#8855ff">⚡ Net Power</text>
  <text x="910" y="506" class="label" fill="#bb88ff">Output</text>
  <line x1="770" y1="490" x2="848" y2="490" stroke="#8855ff" stroke-width="2" stroke-dasharray="6 4"
        style="animation: pipe-flow 1s linear infinite;"/>

  <!-- ════════════════════════════════════════════════════ -->
  <!-- Legend -->
  <!-- ════════════════════════════════════════════════════ -->
  <rect x="900" y="100" width="270" height="180" rx="10" ry="10"
        fill="rgba(10,10,30,0.7)" stroke="rgba(100,120,255,0.15)" stroke-width="1"/>
  <text x="1035" y="125" class="label" fill="#aaa">LEGEND</text>

  <line x1="920" y1="148" x2="960" y2="148" stroke="#8B6914" stroke-width="3" stroke-dasharray="8 5"/>
  <text x="970" y="152" class="label-sm" fill="#ccc" text-anchor="start">Biomass Flow</text>

  <line x1="920" y1="170" x2="960" y2="170" stroke="#00bfff" stroke-width="3" stroke-dasharray="8 5"/>
  <text x="970" y="174" class="label-sm" fill="#ccc" text-anchor="start">Water / Steam (HTC)</text>

  <line x1="920" y1="192" x2="960" y2="192" stroke="#44dd88" stroke-width="3" stroke-dasharray="8 5"/>
  <text x="970" y="196" class="label-sm" fill="#ccc" text-anchor="start">Biogas Flow</text>

  <line x1="920" y1="214" x2="960" y2="214" stroke="#ff6644" stroke-width="3" stroke-dasharray="8 5"/>
  <text x="970" y="218" class="label-sm" fill="#ccc" text-anchor="start">Hot Gas / Combustion</text>

  <line x1="920" y1="236" x2="960" y2="236" stroke="#aabbff" stroke-width="3" stroke-dasharray="8 5"/>
  <text x="970" y="240" class="label-sm" fill="#ccc" text-anchor="start">Air Flow</text>

  <rect x="920" y="254" width="20" height="12" rx="3" stroke="#00d4ff" stroke-width="1.5" fill="rgba(0,212,255,0.15)"/>
  <text x="970" y="264" class="label-sm" fill="#ccc" text-anchor="start">Tank A (Moisture-lean)</text>

  <rect x="920" y="274" width="20" height="12" rx="3" stroke="#00e88f" stroke-width="1.5" fill="rgba(0,232,143,0.15)"/>
  <text x="970" y="284" class="label-sm" fill="#ccc" text-anchor="start">Tank B (Moisture-rich)</text>

  <!-- ── Footer ── -->
  <text x="600" y="625" class="label-sm" fill="#555577">© 2025 Energhx Research Group · Faculty of Engineering, University of Lagos</text>

</svg>
</body>
</html>"""


# ═══════════════════════════════════════════════════════════════════════════════
# CHART BUILDERS  (matplotlib → base64 PNG)
# ═══════════════════════════════════════════════════════════════════════════════

# Dark-theme colours shared by every chart
_BG       = "#0d0d20"
_FACE     = "#0d0d20"
_GRID     = "#1a1a3a"
_TICK     = "#9999bb"
_LABEL    = "#9999bb"
_TITLE    = "#cccccc"
_SPINE    = "#2a2a4a"
_CYAN     = "#00d4ff"
_ORANGE   = "#ff8833"
_TEXT     = "#e0e0f0"


def _fig_to_base64(fig):
    """Render a matplotlib Figure to a data-URI string for html.Img src."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor(), edgecolor="none")
    plt.close(fig)
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("utf-8")
    return f"data:image/png;base64,{b64}"


def _style_ax(ax, xlabel, ylabel):
    """Apply the dark theme to an Axes object."""
    ax.set_facecolor(_FACE)
    ax.set_xlabel(xlabel, color=_LABEL, fontsize=11, fontfamily="sans-serif")
    ax.set_ylabel(ylabel, color=_LABEL, fontsize=11, fontfamily="sans-serif")
    ax.tick_params(colors=_TICK, labelsize=9)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(True, which="major", color=_GRID, linewidth=0.6)
    ax.grid(True, which="minor", color=_GRID, linewidth=0.25, alpha=0.5)
    for spine in ax.spines.values():
        spine.set_color(_SPINE)


def make_hs_chart(states):
    """Create an h-s (enthalpy–entropy) diagram for the HTC steam cycle (matplotlib)."""
    s = states["s"] + [states["s"][0]]
    h = states["h"] + [states["h"][0]]
    labels = states["labels"] + [states["labels"][0]]

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    fig.patch.set_facecolor(_BG)

    # Cycle path
    ax.plot(s, h, color=_CYAN, linewidth=2.4, zorder=2)
    ax.scatter(s, h, s=80, color=_CYAN, edgecolors="#ffffff",
               linewidths=1.5, zorder=3)

    # State-point labels
    for xi, yi, lbl in zip(s, h, labels):
        short = lbl.split("–")[0].strip()
        ax.annotate(short, (xi, yi), textcoords="offset points",
                    xytext=(0, 12), ha="center", fontsize=10,
                    color=_TEXT, fontfamily="monospace", fontweight="bold")

    ax.set_title("h – s Diagram  ·  HTC Steam Cycle",
                 color=_TITLE, fontsize=13, pad=12)
    _style_ax(ax, "Entropy  s  [kJ/(kg·K)]", "Enthalpy  h  [kJ/kg]")
    fig.tight_layout()
    return _fig_to_base64(fig)


def make_t_hdot_chart(states, m_air=1.0):
    """
    Create a T-Ḣ (temperature vs enthalpy-rate) diagram for the gas turbine
    cycle (matplotlib).  Ḣ = ṁ · h  (kW when ṁ in kg/s).
    """
    T = states["T"]
    h = states["h"]
    H_dot = [hi * m_air for hi in h]
    labels = states["labels"]

    # Close the cycle
    T_c = T + [T[0]]
    H_c = H_dot + [H_dot[0]]
    labels_c = labels + [labels[0]]

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    fig.patch.set_facecolor(_BG)

    ax.plot(H_c, T_c, color=_ORANGE, linewidth=2.4, zorder=2)
    ax.scatter(H_c, T_c, s=80, color=_ORANGE, edgecolors="#ffffff",
               linewidths=1.5, zorder=3)

    for xi, yi, lbl in zip(H_c, T_c, labels_c):
        short = lbl.split("–")[0].strip()
        ax.annotate(short, (xi, yi), textcoords="offset points",
                    xytext=(0, 12), ha="center", fontsize=10,
                    color=_TEXT, fontfamily="monospace", fontweight="bold")

    ax.set_title("T – Ḣ Diagram  ·  Gas Turbine Cycle",
                 color=_TITLE, fontsize=13, pad=12)
    _style_ax(ax, "Enthalpy Rate  Ḣ  [kW]", "Temperature  T  [°C]")
    fig.tight_layout()
    return _fig_to_base64(fig)


def empty_chart(title=""):
    """Return a placeholder matplotlib chart as a base64 data-URI."""
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    fig.patch.set_facecolor(_BG)
    ax.set_facecolor(_BG)
    ax.text(0.5, 0.5, "Click  ANALYZE  to generate chart",
            transform=ax.transAxes, ha="center", va="center",
            fontsize=14, color="#555577")
    ax.set_title(title, color="#555577", fontsize=13, pad=12)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.tight_layout()
    return _fig_to_base64(fig)


# ═══════════════════════════════════════════════════════════════════════════════
# DASH APPLICATION
# ═══════════════════════════════════════════════════════════════════════════════

app = dash.Dash(
    __name__,
    title="AD-HTC Gas Cycle Analyzer",
    update_title="Analyzing…",
    meta_tags=[{"name": "viewport",
                "content": "width=device-width, initial-scale=1"}],
)

server = app.server  # for production deployment


# ── Helper to build a parameter slider row ──
def param_row(label, slider_id, min_val, max_val, step, value, unit, marks=None):
    return html.Div(className="param-group", children=[
        html.Div(label, className="param-label"),
        html.Div(f"{value} {unit}", id=f"{slider_id}-display",
                 className="param-value-display"),
        dcc.Slider(
            id=slider_id,
            min=min_val, max=max_val, step=step, value=value,
            marks=marks or {},
            tooltip={"placement": "bottom", "always_visible": False},
        ),
    ])


# ═══════════════════════════════════════════════════════════════════════════════
# LAYOUT
# ═══════════════════════════════════════════════════════════════════════════════

app.layout = html.Div([

    # ── Header ──
    html.Header(className="app-header", children=[
        html.Div([
            html.Div("AD-HTC Fuel-Enhanced Gas Cycle", className="header-title"),
            html.Div("Process Flow Analysis Dashboard", className="header-subtitle"),
        ]),
        html.Div("ENERGHX", className="header-badge"),
    ]),

    # ── Main Content ──
    html.Main(className="main-content", children=[

        # ── 1. SCHEMATIC ──
        html.Div(className="section-header", children=[
            html.Div("🔄", className="section-icon schematic"),
            html.Div([
                html.Div("Process Schematic", className="section-title"),
                html.Div("Animated flow diagram of the AD-HTC fuel-enhanced gas power cycle",
                         className="section-desc"),
            ]),
        ]),
        html.Div(className="glass-card pulse-glow", children=[
            html.Iframe(
                srcDoc=build_schematic_html(),
                className="schematic-frame",
            ),
        ]),

        # ── 2. PARAMETERS ──
        html.Div(className="section-header", children=[
            html.Div("⚙️", className="section-icon params"),
            html.Div([
                html.Div("Input Parameters", className="section-title"),
                html.Div("Configure biomass feed and gas turbine operating conditions",
                         className="section-desc"),
            ]),
        ]),
        html.Div(className="param-grid", children=[

            # Tank A
            html.Div(className="glass-card", children=[
                html.Div("🟦  Tank A — Moisture-Lean Biomass",
                         className="param-card-title tank-a-title"),
                param_row("Feed Rate", "ta-feed", 100, 2000, 50, 500, "kg/h"),
                param_row("Moisture Content", "ta-moist", 5, 50, 1, 20, "%"),
                param_row("HTC Reactor Temp", "ta-treact", 150, 300, 5, 200, "°C"),
            ]),

            # Tank B
            html.Div(className="glass-card", children=[
                html.Div("🟩  Tank B — Moisture-Rich Biomass",
                         className="param-card-title tank-b-title"),
                param_row("Feed Rate", "tb-feed", 100, 2000, 50, 800, "kg/h"),
                param_row("Moisture Content", "tb-moist", 50, 95, 1, 70, "%"),
                param_row("Volatile Solids Fraction", "tb-vs", 0.5, 0.95, 0.05, 0.80, ""),
            ]),

            # Gas Turbine
            html.Div(className="glass-card", children=[
                html.Div("🟧  Gas Turbine Cycle",
                         className="param-card-title turbine-title"),
                param_row("Pressure Ratio", "gt-rp", 4, 25, 1, 10, ""),
                param_row("Turbine Inlet Temp", "gt-tit", 800, 1500, 25, 1200, "°C"),
                param_row("Ambient Temp", "gt-tamb", 10, 45, 1, 25, "°C"),
                param_row("Compressor η_is", "gt-etac", 70, 95, 1, 85, "%"),
                param_row("Turbine η_is", "gt-etat", 70, 95, 1, 90, "%"),
            ]),
        ]),

        # ── Analyze Button ──
        html.Div(className="analyze-btn-wrap", children=[
            html.Button("⚡  Analyze", id="btn-analyze",
                        className="analyze-btn", n_clicks=0),
        ]),

        # ── 3. CHARTS ──
        html.Div(className="section-header", id="charts-section", children=[
            html.Div("📊", className="section-icon charts"),
            html.Div([
                html.Div("Thermodynamic Charts", className="section-title"),
                html.Div("h–s diagram for HTC steam cycle  ·  T–Ḣ diagram for gas turbine cycle",
                         className="section-desc"),
            ]),
        ]),
        html.Div(className="charts-grid", children=[
            html.Div(className="glass-card mpl-chart-wrap", children=[
                html.Img(id="chart-hs-img", className="mpl-chart",
                         src=empty_chart("h – s Diagram")),
            ]),
            html.Div(className="glass-card mpl-chart-wrap", children=[
                html.Img(id="chart-th-img", className="mpl-chart",
                         src=empty_chart("T – Ḣ Diagram")),
            ]),
        ]),

        # ── 4. REPORT ──
        html.Div(className="section-header", children=[
            html.Div("📋", className="section-icon report"),
            html.Div([
                html.Div("Analysis Report", className="section-title"),
                html.Div("State-point data, energy balance, and cycle efficiencies",
                         className="section-desc"),
            ]),
        ]),
        html.Div(id="report-area", className="glass-card", children=[
            html.Div("Configure parameters above and click ANALYZE to generate the report.",
                     style={"textAlign": "center", "padding": "40px",
                            "color": "#555577", "fontSize": "15px"}),
        ]),
    ]),

    # ── Footer ──
    html.Footer(className="app-footer", children=[
        html.Div("© 2025 Energhx Research Group  ·  Faculty of Engineering, University of Lagos",
                 className="footer-text"),
        html.Div("E  N  E  R  G  H  X", className="footer-brand"),
    ]),
])


# ═══════════════════════════════════════════════════════════════════════════════
# CALLBACKS — LIVE SLIDER DISPLAY UPDATES
# ═══════════════════════════════════════════════════════════════════════════════

_slider_units = {
    "ta-feed": "kg/h", "ta-moist": "%", "ta-treact": "°C",
    "tb-feed": "kg/h", "tb-moist": "%", "tb-vs": "",
    "gt-rp": "", "gt-tit": "°C", "gt-tamb": "°C",
    "gt-etac": "%", "gt-etat": "%",
}

for _sid, _unit in _slider_units.items():
    @callback(
        Output(f"{_sid}-display", "children"),
        Input(_sid, "value"),
        prevent_initial_call=False,
    )
    def _update_display(val, unit=_unit):
        return f"{val} {unit}"


# ═══════════════════════════════════════════════════════════════════════════════
# CALLBACK — ANALYZE
# ═══════════════════════════════════════════════════════════════════════════════

@callback(
    Output("chart-hs-img", "src"),
    Output("chart-th-img", "src"),
    Output("report-area", "children"),
    Input("btn-analyze", "n_clicks"),
    State("ta-feed", "value"),  State("ta-moist", "value"),
    State("ta-treact", "value"),
    State("tb-feed", "value"),  State("tb-moist", "value"),
    State("tb-vs", "value"),
    State("gt-rp", "value"),    State("gt-tit", "value"),
    State("gt-tamb", "value"),  State("gt-etac", "value"),
    State("gt-etat", "value"),
    prevent_initial_call=True,
)
def run_analysis(n, ta_feed, ta_moist, ta_treact,
                 tb_feed, tb_moist, tb_vs,
                 gt_rp, gt_tit, gt_tamb, gt_etac, gt_etat):

    # ── Run calculations ──
    gt_states, gt_metrics = brayton_cycle(
        gt_tamb, gt_rp, gt_tit, gt_etac / 100, gt_etat / 100)

    htc_states, htc_metrics = htc_steam_cycle(
        T_reactor=ta_treact, P_reactor=20)

    ad_results = ad_biogas_yield(tb_feed, tb_moist, tb_vs)
    htc_results = htc_process(ta_feed, ta_moist, ta_treact)

    # Estimate air mass flow from biogas energy
    m_air = ad_results["biogas_energy_MJ_h"] / (gt_metrics["q_in"] / 1000) \
        if gt_metrics["q_in"] > 0 else 1.0
    m_air = round(m_air, 2)

    # Net power output
    gt_power_kW = round(gt_metrics["w_net"] * m_air / 3.6, 2)
    htc_power_kW = round(htc_metrics["w_net"] * ta_feed / 3600, 2)
    total_power = round(gt_power_kW + htc_power_kW, 2)

    # ── Build charts ──
    fig_hs = make_hs_chart(htc_states)
    fig_th = make_t_hdot_chart(gt_states, m_air)

    # ── Build report ──
    # Summary metric cards
    metric_cards = html.Div(className="report-grid", children=[
        _metric_card("Net Power", f"{total_power}", "kW"),
        _metric_card("GT Efficiency", f"{gt_metrics['eta_th']}", "%"),
        _metric_card("HTC Efficiency", f"{htc_metrics['eta']}", "%"),
        _metric_card("Biogas Yield", f"{ad_results['biogas_m3_h']}", "m³/h"),
        _metric_card("GT Net Work", f"{gt_metrics['w_net']}", "kJ/kg"),
        _metric_card("Hydrochar", f"{htc_results['hydrochar_kg_h']}", "kg/h"),
        _metric_card("Air Mass Flow", f"{m_air}", "kg/h"),
        _metric_card("Back Work Ratio", f"{gt_metrics['bwr']}", "%"),
    ])

    # Gas Turbine state-point table
    gt_table = _state_table(
        "Gas Turbine Cycle — State Points",
        gt_states, "gt-table",
    )

    # HTC Steam Cycle state-point table
    htc_table = _state_table(
        "HTC Steam Cycle — State Points",
        htc_states, "htc-table",
    )

    # Energy balance table
    energy_data = [
        {"Component": "Compressor Work", "Value": f"{gt_metrics['w_comp']} kJ/kg"},
        {"Component": "Turbine Work", "Value": f"{gt_metrics['w_turb']} kJ/kg"},
        {"Component": "GT Net Work", "Value": f"{gt_metrics['w_net']} kJ/kg"},
        {"Component": "GT Heat Input (Q_in)", "Value": f"{gt_metrics['q_in']} kJ/kg"},
        {"Component": "GT Heat Rejected (Q_out)", "Value": f"{gt_metrics['q_out']} kJ/kg"},
        {"Component": "HTC Steam Pump Work", "Value": f"{htc_metrics['w_pump']} kJ/kg"},
        {"Component": "HTC Steam Turbine Work", "Value": f"{htc_metrics['w_turb']} kJ/kg"},
        {"Component": "HTC Net Work", "Value": f"{htc_metrics['w_net']} kJ/kg"},
        {"Component": "Biogas Energy Output", "Value": f"{ad_results['biogas_energy_MJ_h']} MJ/h"},
        {"Component": "Hydrochar Energy", "Value": f"{htc_results['hydrochar_energy_MJ_h']} MJ/h"},
        {"Component": "HTC Energy Required", "Value": f"{htc_results['energy_required_MJ_h']} MJ/h"},
    ]
    energy_table = html.Div([
        html.Div("Energy Balance", className="report-table-title"),
        dash_table.DataTable(
            columns=[{"name": "Component", "id": "Component"},
                     {"name": "Value", "id": "Value"}],
            data=energy_data,
            style_header={
                "backgroundColor": "rgba(0,212,255,0.08)",
                "color": "#ccc", "fontWeight": "600",
                "border": "1px solid rgba(100,120,255,0.1)",
                "fontFamily": "Inter, sans-serif", "fontSize": "12px",
            },
            style_cell={
                "backgroundColor": "rgba(10,10,30,0.5)",
                "color": "#bbb", "border": "1px solid rgba(100,120,255,0.06)",
                "fontFamily": "JetBrains Mono, monospace", "fontSize": "13px",
                "padding": "10px 14px", "textAlign": "left",
            },
            style_data_conditional=[{
                "if": {"row_index": "odd"},
                "backgroundColor": "rgba(20,20,50,0.5)",
            }],
        ),
    ])

    # AD & HTC process summary
    process_data = [
        {"Parameter": "Tank A Feed Rate", "Value": f"{ta_feed} kg/h"},
        {"Parameter": "Tank A Dry Mass", "Value": f"{htc_results['dry_mass']} kg/h"},
        {"Parameter": "Hydrochar Yield", "Value": f"{htc_results['hydrochar_kg_h']} kg/h"},
        {"Parameter": "Process Water", "Value": f"{htc_results['process_water_kg_h']} kg/h"},
        {"Parameter": "Tank B Feed Rate", "Value": f"{tb_feed} kg/h"},
        {"Parameter": "Tank B Dry Mass", "Value": f"{ad_results['dry_mass']} kg/h"},
        {"Parameter": "Volatile Solids", "Value": f"{ad_results['volatile_solids']} kg/h"},
        {"Parameter": "Biogas Yield", "Value": f"{ad_results['biogas_m3_h']} m³/h"},
        {"Parameter": "Methane Yield", "Value": f"{ad_results['methane_m3_h']} m³/h"},
    ]
    process_table = html.Div([
        html.Div("AD & HTC Process Summary", className="report-table-title"),
        dash_table.DataTable(
            columns=[{"name": "Parameter", "id": "Parameter"},
                     {"name": "Value", "id": "Value"}],
            data=process_data,
            style_header={
                "backgroundColor": "rgba(0,232,143,0.08)",
                "color": "#ccc", "fontWeight": "600",
                "border": "1px solid rgba(100,120,255,0.1)",
                "fontFamily": "Inter, sans-serif", "fontSize": "12px",
            },
            style_cell={
                "backgroundColor": "rgba(10,10,30,0.5)",
                "color": "#bbb", "border": "1px solid rgba(100,120,255,0.06)",
                "fontFamily": "JetBrains Mono, monospace", "fontSize": "13px",
                "padding": "10px 14px", "textAlign": "left",
            },
            style_data_conditional=[{
                "if": {"row_index": "odd"},
                "backgroundColor": "rgba(20,20,50,0.5)",
            }],
        ),
    ])

    report = html.Div([
        metric_cards,
        html.Div(className="report-tables", children=[gt_table, htc_table]),
        html.Br(),
        html.Div(className="report-tables", children=[energy_table, process_table]),
    ])

    return fig_hs, fig_th, report


# ── Report helper functions ──

def _metric_card(label, value, unit=""):
    return html.Div(className="metric-card", children=[
        html.Div(value, className="metric-value"),
        html.Div(label, className="metric-label"),
        html.Div(unit, className="metric-unit"),
    ])


def _state_table(title, states, table_id):
    data = []
    for i in range(len(states["labels"])):
        data.append({
            "State": states["labels"][i],
            "T (°C)": f"{states['T'][i]:.1f}",
            "h (kJ/kg)": f"{states['h'][i]:.1f}",
            "s (kJ/kg·K)": f"{states['s'][i]:.4f}",
        })

    return html.Div([
        html.Div(title, className="report-table-title"),
        dash_table.DataTable(
            id=table_id,
            columns=[
                {"name": "State", "id": "State"},
                {"name": "T (°C)", "id": "T (°C)"},
                {"name": "h (kJ/kg)", "id": "h (kJ/kg)"},
                {"name": "s (kJ/kg·K)", "id": "s (kJ/kg·K)"},
            ],
            data=data,
            style_header={
                "backgroundColor": "rgba(0,212,255,0.08)",
                "color": "#ccc",
                "fontWeight": "600",
                "border": "1px solid rgba(100,120,255,0.1)",
                "fontFamily": "Inter, sans-serif",
                "fontSize": "12px",
            },
            style_cell={
                "backgroundColor": "rgba(10,10,30,0.5)",
                "color": "#bbb",
                "border": "1px solid rgba(100,120,255,0.06)",
                "fontFamily": "JetBrains Mono, monospace",
                "fontSize": "13px",
                "padding": "10px 14px",
                "textAlign": "center",
            },
            style_data_conditional=[
                {"if": {"row_index": "odd"},
                 "backgroundColor": "rgba(20,20,50,0.5)"},
            ],
        ),
    ])


# ═══════════════════════════════════════════════════════════════════════════════
# ENTRY POINT
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("  AD-HTC Fuel-Enhanced Gas Cycle Analyzer")
    print("  Open your browser at  http://127.0.0.1:8050")
    print("=" * 60 + "\n")
    app.run(debug=True, port=8050)
