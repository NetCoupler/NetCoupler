---
title: 'NetCoupler: Inference of Causal Links Between a Network and an External Variable'
tags:
  - R
authors:
  - name: Luke W. Johnston
    orcid: 0000-0003-4169-2616
    equal-contrib: true 
    corresponding: true 
    affiliation: "1, 2"
  - name: Fabian Eichelmann
    orcid: 
    affiliation: 3
  - name: Daniel Ibsen
    orcid: 
    affiliation: 4
  - name: Clemens Wittenbecher
    orcid: 0000-0001-7792-877X
    equal-contrib: true 
    affiliation: 5
affiliations:
 - name: Steno Diabetes Center Aarhus, Denmark
   index: 1
 - name: Aarhus University, Denmark
   index: 2
 - name: 
   index: 3
 - name: 
   index: 4
 - name: SciLifeLab, Division of Food Science and Nutrition, Department of Biology and Biological Engineering, Chalmers University of Technology, Sweden
   index: 5
date: 16 August 2022
bibliography: paper.bib
editor_options: 
  markdown: 
    wrap: 72
---

# Summary

> *A summary describing the high-level functionality and purpose of the
> software for a diverse, non-specialist audience.*

Within the biomedical sciences, data measurement technologies such as
with metabolomics generate large volumes of data. These technologies
usually produce hundreds of variables that analytically can be quite
challenging to adequately analyze. More often, and especially in health
science, researchers are using these data to help answer questions on
what might influence these metabolic characteristics and how these
metabolic characteristics might influence disease states. Essentially,
to begin examining possible causal pathways between "exposures" (like
diet or exercise) and "outcomes" (like diseases such as diabetes).

We created NetCoupler to assist us and our colleagues in answering
research questions that revolve around exploring these possible causal
pathways between an exposure (like diet), a network of hypothetically
connected variables (like proteins within the blood), and an outcome
(like diabetes). While NetCoupler could potentially be used in a wide
variety of research fields and questions, we designed it around largely
datasets that contain many metabolic variables (such as with
metabolomics) and when there is a hypothetical causal pathways
underlying the metabolic variables based on biological rationale.

NetCoupler works by connecting several analytic approaches together in a
way that allows for flexibility in the type of models used throughout
the analysis. For instance, linear regression models can be used when
the data are continuous while survival models can use when it time to
event data. 

# Statement of need

> *Illustrates the research purpose of the software and places it in the
> context of related work*

In fields such as epidemiology, we often work in settings where conducting
controlled experiments to test potential causality is near impossible or at least
extremely difficult. In these areas of research, causal inference or causal
reasoning that doesn't require controlled environments has become more commonly
integrated into study designs and analytic approaches. We also have increasing
access to powerful data generating technologies, such as metabolomics, lipidomics,
or proteomics, and inevitably research questions begin combining causal inference
with these -omic methods. As a result ...

## Underlying principles

## Comparison to other software

# Projects using NetCoupler

The first use of the NetCoupler algorithm (not package) was in Clemens
Wittenbecher's PhD thesis [@Wittenbecher2017]. We have so far presented
work using NetCoupler at several conferences
[@Johnston2020,@Johnston2021a,@Johnston2021,@Johnston2020a] and a
published paper [@Wittenbecher2022]. Aside from these scientific outputs,
NetCoupler is in use in several ongoing projects, such as one in the [UK
Biobank](https://gitlab.com/lwjohnst/ecc-cmd-ukb) and TODO: add Daniel's
project.

# Acknowledgements

Luke Johnston was supported by a Danish Diabetes Academy Postdoctoral
Fellowship Award... TODO: fill in later.

# References
