# Exploring Robustness of ML Fluid Simulation with Algorithmic Improvement

## Introduction

Detailed splashes are a necessity for high-resolution liquid simulations. Modelling large amounts of realistic droplet formation requires the use of very fine spatial discretization due to their complex small scale surface geometry and dynamics. This, in turn, leads to very high computational cost, making it challenging to generate vivid splashes in liquid simulations.

<p align="center">
<img src="https://github.com/vishrutjetly/CSCI596-FinalProject/blob/main/graphics/low_res.gif" width="350" height="350" hspace="50"/>
<img src="https://github.com/vishrutjetly/CSCI596-FinalProject/blob/main/graphics/high_res.gif" width="350" height="350" />
</p>

In the above images, the left has been simulated on a grid size of 64x64, whereas the right has been simulated on a grid size of 512x512. It can be clearly seen that the image on the right is able to render small scale fluid interactions far better compared to the one on left. However, it took only 7 secs for the left image to be simulated, whereas it took roughly 5 mins to simulate the image on the right.

## Methodology

FLIP (Fluid-Implicit-Particle) is a hybrid method of particle and grids, which is currently being widely used in visual effects to model liquid simulations. To improve the visual fidelity of liquid simulations with small-scale details, Kiwon et al. [[1]](#1) utilized neural networks to learn the small-scale splashes from physically accurate simulations. This helped in approximating the sub-grid scale effects that lead to droplet generations, allowing to produce realistic splashes even in coarse simulations and thereby significantly reducing the computational cost. 

In this project, I extend their work by converting their code [[2]](#2) to TensorFlow v2, and studied how neural network learns to model the fluid droplet behaviour with various changes in different parameters. Then I applied Sharpness Aware Minimization to see if it can improve the current model.

### Neural Network Model

The input to the neural network is a feature vector of dimension 27x1, which contains information about the flow. This input is passed through a fully connected layer of size 34, which then produces three outputs - 1) Decision whether this region will turn into a splash or not, 2) Mean of the velocity modification and 3) Variance of the velocity modification.

### Training Data

We know that Neural Networks require a huge amount of training data. In this problem, this is not a huge challenge because we can generate an infinite amount of training data by running fluid simulations over and over with randomized initial positions and velocities. Hence, even changing the parameters within this simulation of fluids for training will also affect the way in which neural networks learn about fluid behaviour. In the following sections, I discuss the effect on the neural network by changing some of these parameters.

### Sharpness Aware Minimization

Machine learning models learn trends in data by minimizing a loss value. With an increase in hyperparameters, the loss function typically gets more complex and non-convex, with multiple local minima and it becomes challenging to achieve generalization. Sharpness Aware Minimization (SAM) is a recently developed procedure that improves model generalization by simultaneously minimizing loss value and loss sharpness. SAM functions by seeking parameters that lie in neighbourhoods having uniformly low loss values. It can be easily implemented [[3]](#3) and I use this method to see if there could be improvements in the fluid simulation.


## Results and Insights

<p align="center">
<img src="https://github.com/vishrutjetly/CSCI596-FinalProject/blob/main/graphics/Flip_vs_MLFlip.gif" style="width:100%">
<font size="24"><strong>Figure 1.</strong> Comparision between FLIP and the Machine Learning FLIP</font>
</p>

Figure 1. shows the comparison between standard flip and machine learning enhanced flip with the same base grid size resolution. Visual improvement can be clearly seen in the case of MLFLIP. Further details are reported in [[1]](#1) 

<br>
<br>

<p align="center">
<img src="https://github.com/vishrutjetly/CSCI596-FinalProject/blob/main/graphics/CompareEpochs.gif" style="width:100%">
<font size="24"><strong>Figure 2.</strong> Improvement in Fluid simulations with increasing epochs</font>

</p>

Figure 2. shows how the fluid simulation improves with increasing epochs. At the start, the fluid behaviour is completely random. This behaviour gradually improves and stabilizes at about 10000 epochs. Beyond this, the visual fidelity almost remains similar and it becomes very difficult to quantitatively measure the accuracy. 

<br>
<br>

<p align="center">
<img src="https://github.com/vishrutjetly/CSCI596-FinalProject/blob/main/graphics/BigData.gif" style="width:100%">
<font size="24"><strong>Figure 3.</strong> Change in neural network's output on increasing the amount of training data by 10 times</font>

</p>

Figure 3. shows the improvement in the simulation by increasing the training data. As stated earlier, this training data can be easily generated by simulating with randomized initial velocities and positions of droplets. Improvement in fluid simulation can be seen in epochs 100 and 1000. Thus, on increasing the training data, the model converges a lot quicker.  

<br>
<br>

<p align="center">
<img src="https://github.com/vishrutjetly/CSCI596-FinalProject/blob/main/graphics/HighResData.gif" style="width:100%">
<font size="24"><strong>Figure 4.</strong> Change in neural network's output on increasing the resolution of training data by a factor of 3</font>

</p>

Figure 4. shows the change in the output when the neural network is trained on simulations with higher resolution (smaller grid size). Almost the same results are obtained in both these cases, hence it can be said that increasing the accuracy of training simulations beyond a certain value does not bring an improvement in the model.

<br>
<br>

<p align="center">
<img src="https://github.com/vishrutjetly/CSCI596-FinalProject/blob/main/graphics/ActualGravity.gif" style="width:100%">
<font size="24"><strong>Figure 5.</strong> Change in neural network's output when trained on simulations in different gravity</font>

</p>

Figure 5. shows the comparison in the output when the neural network is trained on simulations without gravity and with real gravity (9.8). Although it seems training on simulations in zero gravity would result in more splashes, but strangely gravity does not really affect the final outcome.  

<br>
<br>

<p align="center">
<img src="https://github.com/vishrutjetly/CSCI596-FinalProject/blob/main/graphics/NormalVsSAM.gif" style="width:100%">
<font size="24"><strong>Figure 6.</strong> Comparision of normal MLFLIP and MLFLIP improved with SAM</font>

</p>

Figure 6. shows the comparison between the MLFLIP and the improvement obtained on using sharpness aware minimization technique. It can be seen visually how quickly the model converges. MLFLIP with SAM is able to reproduce accurate looking results even at just 100 epochs. Also, the overall visual fidelity seems better when SAM is used.  


## Conclusion and Future Work

The MLFLIP code implemented in [[2]](#2) is converted and updated to TensorFlow v2. Results obtained show in-depth how the machine learning based FLIP is better than standard FLIP simulations. Also, we can visually see how the model evolves and learns with increasing epochs. As expected, the neural network is able to perform better when more training data is fed. However, it was surprising to see little to no improvements in changing the resolution scale and gravity of the training simulations.

This code has been written just to simulate fluids in two dimensions, however, it can be easily updated to incorporate 3D simulations. 
Currently, it is difficult to compare the outputs quantitatively. Hence, we should find some metrics which can provide conclusive confidence in the accuracy of the model.  

## Acknowledgement

I would like to thank Hikaru Ibayashi, Ph.D. student in Computer Science at the University of Southern California, for discussions about the fluid simulations and Sharpness Aware Minimization in the context of this problem, and for his valuable comments after careful reading of this README.


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