# AD-HTC Fuel-Enhanced Gas Cycle Analyzer

An interactive **Dash** web dashboard for analyzing an Anaerobic Digestion and Hydrothermal Carbonization (AD-HTC) fuel-enhanced gas power cycle.

> Developed by **Energhx Research Group** · Faculty of Engineering, University of Lagos

---

## Features

| Feature | Description |
|---------|-------------|
| **Animated SVG Schematic** | Full process-flow diagram with animated pipes, tanks, and rotating turbine |
| **Dual Biomass Tanks** | Tank A (moisture-lean → HTC) and Tank B (moisture-rich → AD) with adjustable parameters |
| **Brayton Cycle Analysis** | Gas turbine thermodynamic calculations with configurable pressure ratio, TIT, and efficiencies |
| **HTC Steam Cycle** | Rankine-like cycle analysis for the HTC reactor subsystem |
| **AD Biogas Yield** | Estimates biogas and methane production from moisture-rich biomass |
| **Matplotlib Charts** | Dark-themed h–s diagram (HTC steam cycle) and T–Ḣ diagram (gas turbine cycle) |
| **Report Generation** | State-point tables, energy balance, and summary metric cards |

---

## Quick Start

### 1. Install Dependencies

```bash
pip install -r requirements.txt
```

### 2. Run the Application

```bash
python app.py
```

### 3. Open in Browser

Navigate to **http://127.0.0.1:8050**

---

## Usage

1. **View the Schematic** — The animated SVG at the top shows the full AD-HTC process flow
2. **Configure Parameters** — Adjust sliders for Tank A, Tank B, and Gas Turbine settings
3. **Click ⚡ Analyze** — Generates thermodynamic charts and a detailed report
4. **Review Results** — Scroll down to see matplotlib charts and the analysis report

---

## Project Structure

```
adhtc-system/
├── app.py              # Main Dash application (layout, callbacks, thermodynamics)
├── requirements.txt    # Python dependencies
├── README.md           # This file
└── assets/
    └── style.css       # Dashboard styling (dark theme)
```

---

## Dependencies

- **dash** ≥ 2.14.0 — Web framework
- **plotly** ≥ 5.18.0 — Required by Dash
- **numpy** ≥ 1.24.0 — Numerical computations
- **matplotlib** ≥ 3.7.0 — Chart rendering (h–s and T–Ḣ diagrams)

---

## Thermodynamic Models

### Brayton (Gas Turbine) Cycle
Calculates compressor/turbine work, heat input/rejection, thermal efficiency, and back-work ratio using isentropic relations with user-defined efficiencies.

### HTC Steam (Rankine-like) Cycle
Simplified Rankine cycle driven by HTC reactor heat, with pump work, turbine work, and cycle efficiency.

### Anaerobic Digestion
Estimates biogas yield from volatile solids content, assuming 60% methane fraction and typical specific yield of 0.40 m³/kg VS.

### HTC Process
Models hydrochar production (60% mass yield) and process energy requirements based on reactor temperature.

---

## Group Members

- Adeyemo Joshua Oluwafemi
- Okoro Joshua Ayomide
- Onwubu Chidubem Rapheal
- Wilcox Samuel
- Adebayo Michael Adewole

---

## License

© 2025 Energhx Research Group · Faculty of Engineering, University of Lagos
