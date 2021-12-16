# CSCI596-FinalProject

## Introduction

Detailed splashes are a necessity for high resolution liquid simulations. Modeling large amounts of realistic droplet formation requires the use of very fine spatial discretization due to their complex small scale surface geometry and dynamics. This in turn, leads to very high computational cost, making it challenging to generate vivid splashes in liquid simulations.


<p float="center">
<img src="https://github.com/vishrutjetly/CSCI596-FinalProject/blob/main/graphics/low_res.gif" width="350" height="350" hspace="50"/>
<img src="https://github.com/vishrutjetly/CSCI596-FinalProject/blob/main/graphics/high_res.gif" width="350" height="350" />
</p>

In the above images, the left has been simulated on a grid size of 64x64, whereas the right has been simulated on a grid size of 512x512. It can be clearly seen that the image on the right is able to render small scale fluid interations far better compared to the one on left. However, it took only 7 secs for the left image to be simulated, whereas it took roughly 5 mins to simulate the image on the right.

## Methodology

FLIP (Fluid-Implicit-Particle) is a hybrid method of particle and grids, which is currently being widely used in visual effects to model liquid simulations. To improve the fisual videlity of liquid simulations with small-scale details, Kiwon *et al* [[1]](#1) utilized neural networks to learn the small-scale splashes from physically accurate simulations. This helped in approximating the sub-grid scale effects that lead to droplet generations, allowing to produce realistic splashes even in coarse simulations and thereby significantly reducing the computational cost. 

In this project, I extend their work by converting their code [[2]](#2) to TensorFlow v2, and studied how neural network learns to model the fluid droplet behaviour with various changes in different parameters. Then I applied Sharpness Aware Minimization to see if it can improve the current model.

### Neural Network Model

The input to the neural network is a feature vector of dimension 27x1, which contains information about the flow. This input is passed through a fully connected layer of size 34, which then produces three outputs - 1) Decision whether this region will turn into splash or not, 2) Mean of the velocity modification and 3) Variance of the velocity modification.

### Training Data

We know that Neural Networks require a huge amount of training data. In this problem, this is not a huge challenge because we can generate infinite amount of training data by running fluid simulations over and over with randomized initial position and velocities. Hence, even changing the parameters within this simulation of fluids for training will also effect the way in which neural networks learn about the fluid behaviour. In following sections, I discuss the effect on the neural network by changing some of these parameters.

### Sharpness Aware Minimization

Machine learning models learn trends in data by minimizing a loss value. With increase in hyperparameters the loss function typically gets more complex and non-convex, with multiple local minima and it becomes challenging to achieve generalization. Sharpness Aware Minimization (SAM) is a recently developed procedure which improves model generalization by simultaneously minimizing loss value and loss sharpness. SAM functions by seeking parameters that lie in neighborhoods having uniformly low loss value. It can be easily implemented [[3]](#3) and I use this method to see if there could be improvements in the fluid simulation.


## Results and Insights

<p float="center">
<figure>
<img src="https://github.com/vishrutjetly/CSCI596-FinalProject/blob/main/graphics/Flip_vs_MLFlip.gif" style="width:100%">
<center><figcaption align = "center">Figure 1. Comparision between FLIP and the Machine Learning FLIP</figcaption></center>
</figure>
</p>

Figure 1. shows the comparision between standard flip and machine learning enhanced flip with same base grid size resolution. Visual improvement can be clearly seen in the case of MLFLIP. Further details are reported in [[1]](#1) 

<br>
<br>
<br>
<br>

<p float="center">
<figure>
<img src="https://github.com/vishrutjetly/CSCI596-FinalProject/blob/main/graphics/CompareEpochs.gif" style="width:100%">
<center><figcaption align = "center">Figure 2. Improvement in Fluid simulations with increasing epochs</figcaption></center>
</figure>
</p>

Figure 2. shows how the fluid simulation improves with increasing epochs. At the start, the fluid behaviour is completely random. This behaviour gradually improves and stabilizes at about 10000 epochs. Beyond this, the visual fidelity almost remains similar and it becomes very difficult to quantitatively measure the accuaracy. 

<br>
<br>
<br>
<br>


<p float="center">
<figure>
<img src="https://github.com/vishrutjetly/CSCI596-FinalProject/blob/main/graphics/BigData.gif" style="width:100%">
<center><figcaption align = "center">Figure 3. Change in neural network's output on increasing the amount of training data by 10 times</figcaption></center>
</figure>
</p>

Figure 3. shows the improvement in the simulation by increasing the training data. As stated earlier, this training data can be easily generated by simulating with randomized initial velocities and positions of droplets. Improvement in fluid simulation can be seen in epochs 100 and 1000. Thus, on increasing the training data, the model converges a lot quicker.  

<br>
<br>
<br>
<br>


<p float="center">
<figure>
<img src="https://github.com/vishrutjetly/CSCI596-FinalProject/blob/main/graphics/HighResData.gif" style="width:100%">
<center><figcaption align = "center">Figure 4. Change in neural network's output on increasing the resolution of training data by a factor of 3</figcaption></center>
</figure>
</p>

Figure 4. shows the change in the output when the neural network is trained on simlations with higher resolution (samller grid size). Almost same results are obtained in both these cases, hence it can be said that increasing the accuracy of training simulations beyond a certain value does not bring an improvement in the model.


<br>
<br>
<br>
<br>


<p float="center">
<figure>
<img src="https://github.com/vishrutjetly/CSCI596-FinalProject/blob/main/graphics/ActualGravity.gif" style="width:100%">
<center><figcaption align = "center">Figure 5. Change in neural network's output when trained on simulations in different gravity</figcaption></center>
</figure>
</p>

Figure 5. shows the comparision in the output when the neural network is trained on simulations without gravity and with real gravity (9.8). Although it seems training on simulations in zero gravity would result in more number of splashes, but strangely gravity does not really effect the final outcome.  

<br>
<br>
<br>
<br>

<p float="center">
<figure>
<img src="https://github.com/vishrutjetly/CSCI596-FinalProject/blob/main/graphics/NormalVsSAM.gif" style="width:100%">
<center><figcaption align = "center">Figure 6. Comparision of normal MLFLIP and MLFLIP improved with SAM</figcaption></center>
</figure>
</p>

Figure 6. shows the comparision between the MLFLIP and the improvement obtained on using sharpness aware minimization technique. It can be seen visually how quickly the model converges. MLFLIP with SAM is able to reproduce accurate looking results even at just 100 epochs. Also, the overall visual fidelity seems better when SAM is used.  


## Conclusion and Future Work

The MLFLIP code implemented in [[2]](#2) is converted and updated to TensorFlow v2. Results obtained show in depth how the machine learning based FLIP is better than standard FLIP simulations. Also, we can visually see how the model evolves and learns with increasing epochs. As expected, neural network is able to perform better when more training data is fed. However it was suprising to see little to no improvements on changing the resolution scale and gravity of the training simuations.

This code has been written just to simulate fluids in two dimension, however it can be easily updated to incorporate 3D simulations. 
Currently, it is difficult to compare the outputs quantitatively. Hence, we should find some metric which can provide conclusive confidence in the accuracy of the model.  

## Acknowledgement

I would like to thank Hikaru Ibayashi, Ph.D. student in Computer Science at the University of Southern California, for discussions about the fluid simulations and Sharpness Aware Minimization in context of this problem.  


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