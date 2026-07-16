****The Silicon Virovore****

An integrated cyber-biotic platform designed to engineer novel peptide inhibitors against Human Endogenous Retrovirus K (HERV-K) envelope proteins and monitor their biosynthesis in real time via closed-loop hardware telemetry.
****🧬 Project Overview****

The Silicon Virovore is a multi-disciplinary, end-to-end therapeutic design and physical monitoring ecosystem. 

The project is split into two primary engines:

**In-Silico Molecular Design Engine**: 

Computationally engineers a custom 21-amino-acid helical peptide inhibitor to target and block the critical binding grooves of the HERV-K surface envelope protein (associated with neurodegenerative diseases).

**Cyber-Biotic Hardware Telemetry Loop:** 

A physical micro-controlled bioreactor simulator that ensures peptide stability during biosynthesis. The system utilizes real-time environmental sensors and closed-loop actuation to prevent thermal denaturation (unfolding) of the engineered proteins.
**
🛠️ System Architecture******

[ NCBI API / Data ] ──> [ Python Bio-Parser ] ──> [ Static C Engine (12GB Heap) ]
                                                            │ (Docking & Design)
                                                            ▼
[ Python Telemetry Dashboard ] <── Serial USB ── [ Elegoo Uno R3 Board (Firmware) ]
             │                                              │ (Closed-Loop)
             ▼                                              ▼
    [ Real-Time Plotting ]                     [ DHT11 Sensor & SG90 Servo Vent ]

1. Computational Software Layer
Data Parsing: Python 3.14 utilizing the NCBI API and high-speed text-parsing libraries to fetch and structure HERV-K target sequences.
Execution Engine: High-performance, memory-optimized C Core utilizing a static 12 GB heap array to run docking, sequence alignments, and structural predictions.
Visualization & Analytics: * PyMOL molecular rendering engines to analyze atomic-level hydrophobic contact maps.
Streamlit local analytics dashboard for interactive sequence mutation tracking.

2. Cyber-Physical Hardware Layer
Microcontroller: Elegoo Uno R3 running low-level C++ firmware.
Sensing: DHT11 digital thermistor tracking environmental temperature inside the synthesis chamber.
Actuation: Active closed-loop physical feedback. When temperature thresholds exceed biological stability limits (≥28 C), the microcontroller:

Triggers an SG90 micro-servo to actuate a physical cooling vent.
Illuminates a physical Red LED warning indicator.
Telemetry Link: A Python-based Serial (pyserial) bridge that streams raw physical telemetry out of the microcontroller's USB port directly into a terminal-based live CLI dashboard.

**🚀 Repository Structure**

Code snippet
├── src/
│   ├── parser.py              # NCBI sequence fetcher & data parser
│   ├── main.c                 # C-based execution core
│   └── firmware/
│       └── bioreactor_loop.ino # Elegoo Uno R3 Arduino firmware
├── telemetry/
│   └── virovore_monitor.py    # Python Serial telemetry dashboard
├── models/
│   └── docking.res.1.pdb      # Target receptor & peptide complex coordinates
├── output/
│   └── binding_pose_1.png     # High-resolution PyMOL molecular render
└── app.py                     # Streamlit local analytics dashboard

**📈 Real-Time Telemetry Interface**

When the physical hardware loop is active, the Python telemetry engine outputs live structural integrity metrics:

Bash
==============================================================
VIROVORE BIOREACTOR CYBER-BIOTIC MONITORING SYSTEM
==============================================================
[*] Successfully connected to hardware on port: /dev/cu.usbmodem14101
[*] Waiting for incoming bioreactor telemetry...

[🟢 SAFE ] Temp: 24.20°C | Bioreactor Vent: CLOSED [Stable]
[🟢 SAFE ] Temp: 24.50°C | Bioreactor Vent: CLOSED [Stable]
[🚨 ALERT] Temp: 28.10°C | Bioreactor Vent: OPEN [Venting Heat]
[🚨 ALERT] Temp: 28.40°C | Bioreactor Vent: OPEN [Venting Heat]
🏆 GSEF / International Science and Engineering Fair (ISEF) Compliance
This project is fully compliant with all ISEF Rules and Regulations:
100% In-Silico & Mock-Hardware Design: No physical viral agents, live pathogens, or biological materials are cultured or handled. The project bypasses hazardous biological agent restrictions (no BSL-2 wet-lab approval required).
Low-Voltage Electrical Safety: The physical monitoring loop is completely powered by low-voltage 5V DC via USB, ensuring 100% compliance with display safety rules.
