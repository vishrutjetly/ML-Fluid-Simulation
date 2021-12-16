# CSCI596-FinalProject

## Introduction

Detailed splashes are a necessity for high resolution liquid simulations. Modeling large amounts of realistic droplet formation requires the use of very fine spatial discretization due to their complex small scale surface geometry and dynamics. This in turn, leads to very high computational cost, making it challenging to generate vivid splashes in liquid simulations.


<p float="center">
<img src="https://github.com/vishrutjetly/CSCI596-FinalProject/blob/main/graphics/low_res.gif" width="350" height="350" hspace="50"/>
<img src="https://github.com/vishrutjetly/CSCI596-FinalProject/blob/main/graphics/high_res.gif" width="350" height="350" />
</p>


## Methodology

FLIP (Fluid-Implicit-Particle) is a hybrid method of particle and grids, which is currently being widely used in visual effects to model liquid simulations. To improve the fisual videlity of liquid simulations with small-scale details, Kiwon *et al* [[1]](#1) utilized neural networks to learn the small-scale splashes from physically accurate simulations. This helped in approximating the sub-grid scale effects that lead to droplet generations, allowing to produce realistic splashes even in coarse simulations and thereby significantly reducing the computational cost. 

In this project, I extend their work by converting their code [[2]](#2) to TensorFlow v2, and studied how neural network learns to 

### Neural Network Model

### Training Data

### Sharpness Aware Minimization


## Results and Insights


## Conclusion and Future Work


## Acknowledgement



## References

<ol>
	<li> 
		<a id="1">Kiwon Um</a>, Xiangyu Hu, and Nils Thuerey. "Liquid splash modeling with neural networks." <i>Computer Graphics Forum</i>. Vol. 37. No. 8. 2018.
	</li>
	<li>
		<a id="2">"MLFLIP"</a> (available at: <a href="https://github.com/kiwonum/mlflip">https://github.com/kiwonum/mlflip/</a>) (Accessed: 23 November 2021)
	</li>
	<li>
		<a id="3">Pierre Foret,</a> et al. "Sharpness-aware minimization for efficiently improving generalization." arXiv preprint arXiv:2010.01412 (2020).
	</li>
	<li>
		<a id="4">"Mantaflow</a>-An extensible framework for fluid simulation" (available at: <a href="http://mantaflow.com/">http://mantaflow.com/</a>) (Accessed: 23 November 2021)
	</li>
</ol>