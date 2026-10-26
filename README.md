# Inversión Acústica de Onda Completa mediante Redes Neuronales Informadas por la Ecuación Eikonal sobre Señales Crudas de Radiofrecuencia Ultrasónica

[![PyTorch](https://img.shields.io/badge/PyTorch-2.x%2B-EE4C2C?logo=pytorch&logoColor=white)](https://pytorch.org/)
[![CUDA](https://img.shields.io/badge/CUDA-Accelerated-76B900?logo=nvidia&logoColor=white)](https://developer.nvidia.com/cuda-zone)
[![FastAPI](https://img.shields.io/badge/FastAPI-0.115%2B-009688?logo=fastapi&logoColor=white)](https://fastapi.tiangolo.com/)
[![React](https://img.shields.io/badge/React_18-Three.js_%2F_WebGL-61DAFB?logo=react&logoColor=black)](https://react.dev/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

> **Autor:** Yerson Hernando León-Montañez ([ORCID: 0009-0000-4439-8661](https://orcid.org/0009-0000-4439-8661))  
> *Facultad de Ingeniería, Maestría en Inteligencia Artificial, Universidad de La Salle, Bogotá, Colombia.*

---

## 1. Descripción Física y Matemática del Marco PINN-Eikonal

La reconstrucción cuantitativa simultánea de la **velocidad de fase compresional** $c(\mathbf{x})$ y del **coeficiente de atenuación espectral** $\alpha_0(\mathbf{x})$ en medios viscoelásticos heterogéneos a partir de mediciones ultrasónicas en reflexión unilateral (*pulse-echo*) es un problema inverso no lineal severamente mal planteado debido a la abertura angular limitada y al acoplamiento cruzado (*cross-talk*) con el patrón de interferencia coherente (*speckle*).

Este repositorio implementa un marco variacional libre de malla basado en **Redes Neuronales de Representación Implícita Periódica (SIREN)** informadas por la **ecuación diferencial parcial Eikonal factorizada** y el **transporte log-espectral regularizado**, acopladas directamente a operadores físicos diferenciables sobre señales crudas de radiofrecuencia (RF):

$$\lvert\nabla\tau(\mathbf{x})\rvert^2 = \left(\frac{\partial\tau}{\partial x}\right)^2 + \left(\frac{\partial\tau}{\partial z}\right)^2 = \frac{1}{c^2(\mathbf{x})}, \qquad 2\,\nabla\tau(\mathbf{x})\cdot\nabla\ln A(\mathbf{x}, f) + \nabla^2\tau(\mathbf{x}) = -\frac{2\,\alpha_0(\mathbf{x})\,f}{c(\mathbf{x})}$$

### Componentes Principales de la Arquitectura

1. **`DecoupledUltrasoundPINN` ([`models/pinn_net.py`](models/pinn_net.py)):**
   - Dos sub-redes SIREN estrictamente independientes ($\nabla_\theta \alpha_\psi \equiv \mathbf{0}$ y $\nabla_\psi c_\theta \equiv \mathbf{0}$) con frecuencias portadoras diferenciadas:
     - **Rama Cinemática ($c_\theta$):** $\omega_{0,c} = 6.5$, proyectada al intervalo físico $[1460, 1650]\text{ m/s}$.
     - **Rama Disipativa ($\alpha_\psi$):** $\omega_{0,\alpha} = 7.8$, proyectada al intervalo $[0.0, 1.35]\text{ dB/MHz/cm}$.
2. **`FactorizedEikonalTravelTimeNet` ([`physics/diff_das_eikonal.py`](physics/diff_das_eikonal.py)):**
   - Descompone el tiempo de viaje como $\tau_\phi(\mathbf{x}) = \tau_0(\mathbf{x}) + \mathcal{W}(z)\,\delta\tau_\phi(\tilde{\mathbf{x}})$ con $\tau_0(\mathbf{x}) = z/c_0$, eliminando analíticamente la singularidad de fuente en el plano del transductor $z=0$.
3. **`DifferentiablePlaneWaveDAS` — Régimen I (Ondas Planas Multi-Ángulo Pre-Formadas, CUBDL):**
   - Operador Delay-and-Sum (DAS) analítico y diferenciable en GPU sobre tensores RF de $N_{\text{ang}} = 75$ ángulos ($\theta_n \in [-16^\circ, +16^\circ]$) y $N_{\text{el}} = 128$ canales con descomposición en sub-aperturas angulares ($\Theta_L < -2^\circ$ vs. $\Theta_R > +2^\circ$) y maximización del Factor de Coherencia de Mallart-Fink ($\text{CF}$).
4. **Operador de Pantalla de Fase en Número de Onda Lateral ($k_x$) y Derivada Axial en Escala de Fresnel — Régimen II (Campos Focalizados Post-Formados, OASBUD):**
   - Descomposición sintética mediante FFT lateral $\mathcal{F}_x\{s_a(z,x)\}$ en sub-aperturas opuestas ($k_x < 0$ vs. $k_x > 0$), re-enfoque diferenciable de pantalla de fase (*Phase Screen*), normalización lateral por cuantil de profundidad $Q_{0.78}(z)$ para suprimir capas especulares horizontales (epidermis y ligamentos de Cooper), derivada axial convolucional Gaussiana a escala de la primera zona de Fresnel y calibración viscoelástica causal mediante las relaciones de **Beer-Lambert y Kramers-Kronig**.

---

## 2. Resumen de Métricas Cuantitativas Validadas

### Fase 1: Benchmark Internacional CUBDL (Ondas Planas Multi-Ángulo Pre-Formadas, $N=7$)

| Dataset | Medio Experimental | Centro $(x_t, z_t)$ [mm] | $\text{CNR}_c$ [dBr] | $\text{CNR}_\alpha$ [dBr] | $d_{\text{pico},c}$ [mm] | $d_{\text{pico},\alpha}$ [mm] | $\Delta\text{CF}$ [%] | Residuo Eikonal $\mathcal{L}_{\text{Eikonal}}$ |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`MYO001`** | Fantoma CIRS 040GSE | $(+2.0,\, 18.5)$ | $9.82$ | $9.91$ | $1.33$ | $0.75$ | $+20.30\%$ | $3.15 \times 10^{-5}$ |
| **`MYO002`** | Fantoma CIRS 040GSE | $(+0.0,\, 23.5)$ | $9.40$ | $9.38$ | $0.76$ | $0.25$ | $+14.18\%$ | $2.21 \times 10^{-5}$ |
| **`MYO003`** | Fantoma CIRS 040GSE (borde) | $(+16.0,\, 18.0)$ | $9.21$ | $7.66$ | $2.53$ | $2.41$ | $+19.86\%$ | $2.41 \times 10^{-5}$ |
| **`MYO004`** | Fantoma CIRS 040GSE | $(-2.5,\, 11.0)$ | $9.88$ | $9.66$ | $1.62$ | $1.03$ | $+34.91\%$ | $3.24 \times 10^{-5}$ |
| **`MYO005`** | Fantoma CIRS 040GSE | $(+4.0,\, 11.5)$ | $9.46$ | $9.91$ | $1.56$ | $1.13$ | $+35.29\%$ | $3.13 \times 10^{-5}$ |
| **`MYO006`** | Fantoma CIRS 040GSE (borde) | $(-16.5,\, 10.5)$ | $8.04$ | $8.70$ | $2.07$ | $2.36$ | $+29.01\%$ | $2.76 \times 10^{-5}$ |
| **`TSH002`** | Tejido mamario *in-vivo* | $(-6.5,\, 20.0)$ | $10.12$ | $10.20$ | $0.50$ | $0.31$ | $+4.61\%$ | $3.17 \times 10^{-5}$ |
| **Media $\pm$ DE** | **Serie CIRS (`MYO001`–`006`)** | **Radio $R_0 = 3.0\text{ mm}$** | **$9.30 \pm 0.61$** | **$9.20 \pm 0.80$** | **$1.65 \pm 0.56$** | **$1.32 \pm 0.81$** | **$+25.59\%$** | **$2.82 \times 10^{-5}$** |

### Fase 2: Repositorio Clínico OASBUD (Inversión de Campo Completo 100% Ciega, $f_s = 40\text{ MHz}$)

| Caso RF (OASBUD) | Régimen Físico | $\bar{c}_{\text{ROI}} \,/\, c_{\max}$ [m/s] | $\bar{c}_{\text{bg}} \pm \sigma_c$ [m/s] | $\bar{\alpha}_{0,\text{ROI}} \,/\, \alpha_{\max}$ [dB/MHz/cm] | $\text{CNR}_{c \,/\, \alpha}$ [dBr] | Dice Ciego $(c \,/\, \alpha_0)$ | IoU Ciego $(c \,/\, \alpha_0)$ | $d_{\text{pico}, c \,/\, \alpha}$ [mm] | $\Delta\text{CF}_{k_x}$ [%] |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`OASBUD_1mm`** | Baja dispersión (BI-RADS 3) | $1549.7 \,/\, 1554.0$ | $1541.3 \pm 2.2$ | $0.59 \,/\, 0.64$ | $14.26 \,/\, 16.87$ | $0.760 \,/\, 0.867$ | $0.612 \,/\, 0.766$ | $1.05 \,/\, 1.02$ | $+23.98\%$ |
| **`OASBUD_2mm`** | Baja dispersión (BI-RADS 3) | $1552.1 \,/\, 1557.3$ | $1541.2 \pm 2.4$ | $0.62 \,/\, 0.66$ | $15.46 \,/\, 17.40$ | $0.827 \,/\, 0.859$ | $0.705 \,/\, 0.752$ | $0.97 \,/\, 0.64$ | $+28.61\%$ |
| **`OASBUD_4lw`** | Alta dispersión (BI-RADS 5) | $1565.9 \,/\, 1572.8$ | $1543.8 \pm 5.7$ | $0.76 \,/\, 0.86$ | $13.48 \,/\, 14.67$ | $0.787 \,/\, 0.748$ | $0.649 \,/\, 0.598$ | $1.89 \,/\, 3.53$ | $+36.87\%$ |
| **`OASBUD_8he`** | Alta dispersión (BI-RADS 5) | $1568.0 \,/\, 1581.3$ | $1544.8 \pm 7.7$ | $0.80 \,/\, 0.96$ | $12.50 \,/\, 14.14$ | $0.713 \,/\, 0.709$ | $0.554 \,/\, 0.550$ | $1.89 \,/\, 1.98$ | $+42.62\%$ |
| **Promedio** | **Baja dispersión (Clase 0)** | **$1550.9 \,/\, 1555.7$** | **$1541.3 \pm 2.3$** | **$0.61 \,/\, 0.65$** | **$14.86 \,/\, 17.14$** | **$0.794 \,/\, 0.863$** | **$0.659 \,/\, 0.759$** | **$1.01 \,/\, 0.83$** | **$+26.30\%$** |
| **Promedio** | **Alta dispersión (Clase 1)** | **$1567.0 \,/\, 1577.1$** | **$1544.3 \pm 6.7$** | **$0.78 \,/\, 0.91$** | **$12.99 \,/\, 14.41$** | **$0.750 \,/\, 0.729$** | **$0.602 \,/\, 0.574$** | **$1.89 \,/\, 2.76$** | **$+39.75\%$** |

---

## 3. Estructura del Repositorio

```text
pinn_ultrasonido_helmholtz/
├── main_eikonal_pinn_cubdl.py       # Motor de inversión PINN-Eikonal + DAS diferenciable (Régimen I: CUBDL)
├── main_eikonal_pinn_oasbud.py      # Motor de inversión PINN-Eikonal en campo completo (Régimen II: OASBUD)
├── models/
│   └── pinn_net.py                  # Arquitectura implícita SIREN desacoplada (DecoupledUltrasoundPINN)
├── physics/
│   ├── diff_das_eikonal.py          # Red Eikonal factorizada y operador DifferentiablePlaneWaveDAS
│   └── helmholtz_loss.py            # Operadores diferenciales y funcionales variacionales
├── data_engine/
│   └── loader.py                    # Cargadores de tensores RF crudos HDF5 (CUBDL) y MAT (OASBUD)
├── medical_insight_webapp/          # Plataforma web interactiva Medical Insight 3D
│   ├── backend/                     # Servidor API REST y gestor de trabajos GPU en FastAPI
│   └── frontend/                    # Cliente React 18 + TypeScript + Three.js + TailwindCSS
├── start_webapp.ps1 / .bat          # Scripts de arranque unificado para Windows
└── requirements.txt                 # Dependencias científicas de Python y PyTorch
```

---

## 4. Instalación, Descarga de Datasets y Ejecución

### 5.1. Requisitos e Instalación del Entorno Python

Se recomienda Python 3.10+ y una GPU compatible con NVIDIA CUDA:

```bash
git clone https://github.com/yhleonm/pinn_ultrasonido_helmholtz.git
cd pinn_ultrasonido_helmholtz
python -m venv .venv
# En Windows (PowerShell):
.\.venv\Scripts\Activate.ps1
# En Linux / macOS:
# source .venv/bin/activate

pip install -r requirements.txt
```

### 5.2. Descarga de Datasets Públicos de Radiofrecuencia (RF)

Por buenas prácticas de ingeniería de repositorios científicos, los archivos binarios crudos pesados (`*.hdf5`, `*.mat`) están excluidos en `.gitignore`. Descárguelos directamente desde sus fuentes académicas oficiales:

1. **Benchmark CUBDL (Ondas Planas Pre-Formadas, Task 1 `.hdf5`):**
   - Repositorio oficial IEEE DataPort / GitLab CUBDL: [https://gitlab.com/dongwoon.hyun/cubdl](https://gitlab.com/dongwoon.hyun/cubdl) (Referencia: Hyun et al., *IEEE TUFFC*, 2021).
   - Coloque los archivos `MYO001.hdf5`–`MYO006.hdf5` y `TSH002.hdf5` en `data_cubdl_prueba_1/datasets_descargados/`.
2. **Base de Datos Clínica OASBUD (`OASBUD.mat`):**
   - Repositorio oficial Zenodo (Registro **545928**): [https://zenodo.org/records/545928](https://zenodo.org/records/545928) (Referencia: Piotrzkowska-Wróblewska et al., *Medical Physics*, 2017).
   - Descargue `OASBUD.mat` y ubíquelo en `data_raw/OASBUD.mat`.

### 5.3. Ejecución de la Inversión PINN-Eikonal en GPU

- **Para ejecutar el Régimen I (CUBDL — Ondas Planas Multi-Ángulo):**
  ```bash
  python main_eikonal_pinn_cubdl.py --dataset MYO001 --epochs 140
  ```
- **Para ejecutar el Régimen II (OASBUD — Inversión Ciega de Campo Completo a 40 MHz):**
  ```bash
  python main_eikonal_pinn_oasbud.py
  ```

### 4.4. Despliegue del Visor Interactivo *Medical Insight 3D*

Para levantar la plataforma interactiva *Medical Insight 3D*:

```powershell
# Opción rápida en Windows (levanta Backend FastAPI en puerto 8000 y Frontend Vite en puerto 5173):
.\start_webapp.ps1
```

O manualmente en dos terminales:
```bash
# Terminal 1 — Backend FastAPI:
uvicorn medical_insight_webapp.backend.app.main:app --reload --port 8000

# Terminal 2 — Frontend React + Three.js:
cd medical_insight_webapp/frontend
npm install
npm run dev
```

Abra [http://localhost:5173](http://localhost:5173) en su navegador para inspeccionar en 2D/3D todos los casos de **Fase 1 (CUBDL)** y **Fase 2 (OASBUD)**.

---

## 6. Citación (BibTeX)

Si utiliza este código, los operadores diferenciables Eikonal-DAS / Phase Screen $k_x$ o los resultados experimentales en su investigación, por favor cite el artículo:

```bibtex
@article{LeonMontanez2026PINNEikonal,
  author  = {León-Montañez, Yerson Hernando},
  title   = {Inversión acústica de onda completa mediante redes neuronales informadas por la ecuación {Eikonal} sobre señales crudas de radiofrecuencia ultrasónica},
  journal = {Revista de la Academia Colombiana de Ciencias Exactas, Físicas y Naturales},
  year    = {2026},
  note    = {Sección: Ciencias Físicas},
  url     = {https://github.com/yhleonm/pinn_ultrasonido_helmholtz}
}
```
