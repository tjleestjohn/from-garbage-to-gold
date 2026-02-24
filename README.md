# From Garbage to Gold: A Data-Architectural Theory of Predictive Robustness

**Links**

[![Read Paper (Status: Under Review)](https://img.shields.io/badge/Read%20Paper-(Status:%20Under%20Review)-blue?style=for-the-badge&logo=adobeacrobatreader&logoColor=white)](https://arxiv.org/your-paper-link)
[![Paper License](https://img.shields.io/badge/Paper_License-CC_BY_4.0-green?style=for-the-badge)](https://creativecommons.org/licenses/by/4.0/)
[![Code License](https://img.shields.io/badge/Code_License-CC_BY--NC_4.0-lightgrey?style=for-the-badge)](https://creativecommons.org/licenses/by-nc/4.0/)

<p align="center">
  <a href="#usage"><strong>Jump to R Simulation Script 💻</strong></a> &nbsp;|&nbsp; 
  <a href="./G2G_Paper_Audio_Overview.mp3?raw=true"><strong>🎧 Download Audio Overview (MP3)</strong></a>
</p>

### 📰 News
* **Feb 24, 2026:** Paper submitted to arXiv (`cs.LG`).
* **Feb 24, 2026:** Initial code release (v1.0).

> **TL;DR**
>
> We prove mathematically why adding **MORE** error-prone variables can beat cleaning **FEWER** ones to perfection.
>
> *The key*: partitioning predictor-space noise into "Predictor Error" (fixable by cleaning or by adding diverse measurements) and "Structural Uncertainty" (ambiguities due to stochastic generative mappings — only fixable by adding diverse measurements).
>
> This manuscript helps explain why modern ML models can succeed despite being trained on dirty data.

---

## 🏗️ Why This Matters
The AI community is obsessed with foundation models for unstructured data (text, images), but **80%+ of enterprise data is tabular**, stored in massive "data swamps" — vast data warehouses that conventional wisdom says are unusable without cleaning.

Our framework suggests these swamps, if approached correctly, may be untapped gold mines that can be leveraged without having to invest heavily in cleaning pipelines.

## 🎯 Who Should Read This
* **ML practitioners** working with tabular data in any enterprise settings where data cleaning is impractical at scale.
* **Researchers** interested in theories of high-dimensional statistics, benign overfitting, information theory, or data-centric AI.
* **Anyone building ML systems** on real-world, high-dimensional tabular data who's tired of being told to "clean it all first."

## ⚠️ What We Don't Claim
This isn't universal. The framework requires data with a **latent hierarchical structure** (e.g., medical diagnoses driven by unmeasured physiology, stock prices driven by hidden sentiment, sensor readings driven by underlying physical states).

Additionally, this is not a magic pill. Pre-processing effort shifts from ensuring data cleanliness to ensuring robust data architecture. Fortunately, this architectural shift scales significantly better than cleaning in high-dimensional environments.

We detail the boundaries in the paper.

## 💡 Context
As the first author of the paper/simulation, I've spent the past 2.5+ years formalizing this theory to explain anomalies my coauthors and I kept encountering in industry: models trained on error-prone, high-dimensional tabular data sometimes achieve state-of-the-art predictive performance, defying the "Garbage In, Garbage Out" mantra.

We are now operationalizing this theory into a **warehouse-native infrastructure** (e.g., Snowflake, Databricks, etc.) to allow companies to create bespoke, reliable predictive models directly from their enterprise data swamps, saving on cleaning costs (which can consume 60-80% of ML project budgets) while also unlocking value in a previously underutilized resource.

---

## 🧠 Core Thesis (contrarian, yet simple)

1. **Data cleaning has a theoretical limit.** You can clean a fixed set of variables until it is perfect, but if the set collectively only serves as an incomplete *proxy* for the underlying latent drivers of the system, you are still stuck with **"Structural Uncertainty."**
2. **The Fix:** We show that **expanding the set strategically** — even when doing so adds variables that contain errors — allows for efficient, perfect triangulation of the latent drivers. By treating high-dimensional redundancy as a feature rather than a bug, we can actually lower the error rate significantly below that of a smaller set of perfectly clean variables.

## 📉 Consequences

* **Data Architecture over Hygiene:** We argue for a rebalancing of priorities. Instead of obsessing over data cleaning, more emphasis should be put towards building an architecturally sound portfolio of variables, one that contains complete and redundant coverage of the latent drivers.
* **New Feature Selection Strategies:** We propose updated feature selection strategies for big data applications. These methods identify sets of features that efficiently maximize robustness against predictor-space noise, allowing practitioners to remain inside their computational budgets.
* **Principled Support for "Dirty Inputs, Clean Labels":** We offer a data architecture explanation of where and why label cleaning is still distinctly powerful, concluding that allocating cleaning budget almost exclusively to labels, while relying on data architecture to handle predictor-space noise, is a pragmatic strategy when systematic error jointly affects both.
* **Democratizing Deployment (Local Factories over Universal Models):** Rather than deploying pre-trained foundation models that often degrade across sites, our framework supports building "local factories" that learn from site-specific data, turning institutional data swamps into predictive assets.

---

## 📊 Empirical Evidence

### 1. Simulation
This Repo includes a fully-annotated simulation in R that demonstrates the core mechanisms: how "Dirty Breadth" can beat "Clean Parsimony," and how strategic feature selection speeds up convergence.

**[📂 View g2g_simulation.R](./g2g_simulation.R)**

### 2. Real-World Motivation
We discuss the motivating case study with **Cleveland Clinic Abu Dhabi** data where we achieved a **.909 AUC** using 558k+ patients, 3.4+ million patient-months, and thousands of uncurated, error-prone variables, beating the standard-of-care curated models.

**[Read the Clinical Case Study from Cleveland Clinic Abu Dhabi (PLOS Digital Health)](https://journals.plos.org/digitalhealth/article?id=10.1371/journal.pdig.0000589)**

---

## 🔗 A Multidisciplinary Synthesis
Dismantling the universality of the "Garbage In, Garbage Out" dogma requires thorough coverage of many core and adjacent topics. As a result, the paper is long and dense (**110+ pages of text and equations**) with multiple connected threads.

We didn't invent new math or logic though; we synthesized the old with the new.

The paper bridges the gap between modern ML, High-Dimensional Statistics, Benign Overfitting, Psychometrics, Latent Factor Models, Information Theory, and practical Data Engineering **to explain why it is possible to get gold predictions from data traditionally considered garbage.**

---

## 💬 Open Forum
We would love to hear your thoughts on the theory and its potential applications. Happy to answer any questions.

---

## <a id="usage"></a>💻 Usage (R Script)
To downlaod the simulation script, simply copy, paste, and run the below code in R. Or download it from this Repo directly.

When the g2g_simulation.R script opens in R, select all, then run.

Follow the prompts in the R console to proceed through the simulation stages.

**Notes:**
* **Automatic Setup:** The script is self-contained. It will automatically detect and install any missing dependencies (e.g., `data.table`, `ranger`, `ggplot2`) upon the first run.
* **Runtime:** The simulation works in a high-dimensional predictor space with a large sample size. Please expect substantial execution time.

```R
# 1. Define the direct URL to the raw script on GitHub
url <- "https://raw.githubusercontent.com/tjleestjohn/from-garbage-to-gold/main/g2g_simulation.R"

# 2. Define what you want to name the file on your computer
file_name <- "g2g_simulation.R"

# 3. Download just the script
download.file(url, destfile = file_name)

# 4. Open the script in your editor (like RStudio)
file.edit(file_name)
```

---

## 📬 Citation & Contact

If you use this work, please cite the paper:

```bibtex
@article{lee-stjohn2026garbage,
  title={From Garbage to Gold: A Data-Architectural Theory of Predictive Robustness},
  author={Lee-St. John, Terrence J. and Lawson, Jordan L. and Piechowski-Jozwiak, Bartlomiej},
  journal={arXiv preprint arXiv:25XX.XXXXX},
  year={2026},
  url={https://arxiv.org/abs/25XX.XXXXX}
}

```

| **Resource** | **Link** |
| --- | --- |
| **From Garbage to Gold - Full Paper** | [arXiv cs.LG](https://www.google.com/search?q=https://arxiv.org/abs/25XX.XXXXX) |
| **Cleveland Clinic Abu Dhabi Case Study** | [PLOS Digital Health](https://journals.plos.org/digitalhealth/article?id=10.1371/journal.pdig.0000589) |
| **Contact First Author** | [Email Me](mailto:tjleestjohn@gmail.com) |

```
