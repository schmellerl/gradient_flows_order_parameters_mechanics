
<h1>Gradient flows for coupling order parameters and mechanics</h1>

FENICS implementation of the examples from

Schmeller, L. & Peschka, D. (2022). [10.20347/WIAS.PREPRINT.2909](http://dx.doi.org/10.20347/WIAS.PREPRINT.2909)

supplementary simulation data available at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5832662.svg)](https://doi.org/10.5281/zenodo.5832662)

**Abstract:** We construct the formal gradient flow structure for phase field evolution coupled to mechanics in Lagrangian coordinates, present common ways to couple the evolution and discuss general incremental minimization strategies. While the usual presentation of continuum mechanics is intentionally very brief, the focus of this paper is on an extensible functional analytical framework and a discretization approach that preserves an appropriate variational structure as much as possible. As examples, we first present phase separation and swelling of gels and then the approach of stationary states of multiphase systems with surface tension and show the robustness of the general approach.

<h3>1. Introductory Example</h3>

In this introductory example (Section 4.1 in the paper) we examine general aspects of gradient flow discretizations. More details concerning the abstract discretization strategy of gradient flows can be found in the paper. Documentation on the discretization using the *Unified Form Language* (UFL) and FEniCS can be found here: [fenicsproject.org](https://fenicsproject.org).

View documented notebook [Example_Intro.ipynb](colab/Example_Intro.ipynb) or open in Google Colab. [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/schmellerl/gradient_flows_order_parameters_mechanics/blob/main/colab/Example_Intro.ipynb)

<h3>2. Phase separation and gels with Darcy dissipation </h3>

In this example (Section 4.2 in the paper) we present Phase separation using a Flory-Huggins potential on a circular domain. We introduce constraints to the problem using a Lagrange multiplier.
  
View documented notebook [Example2.ipynb](colab/Example2.ipynb) or open in Google Colab. [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/schmellerl/gradient_flows_order_parameters_mechanics/blob/main/colab/Example2.ipynb)

<h3>3. Multiphase systems with Stokes dissipation </h3>

In this example (Section 4.3 in the paper) we show a three phase system using two phase field parameters. 

View documented notebook [Three_Phase_Stokes.ipynb](colab/Three_Phase_Stokes.ipynb) or open in Google Colab. [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/schmellerl/gradient_flows_order_parameters_mechanics/blob/main/colab/Three_Phase_Stokes.ipynb)

View documented notebook [Three_Phase_pureDef.ipynb](colab/Three_Phase_pureDef.ipynb) or open in Google Colab. [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/schmellerl/gradient_flows_order_parameters_mechanics/blob/main/colab/Three_Phase_pureDef.ipynb)
