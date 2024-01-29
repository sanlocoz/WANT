// The hints are written in a dictionary, the key is the question number. In de Rmd thesse numbers appear in the button 'type'
//numbers refer to: first digit->assignment number; second digit->version of hint (started with 0); third and fourth:hint number
//commented dict[] numbers below, give the original number 11/12/2019
dict = {};
//dict[1] = "It can be handy to define the fluxes and boundary conditions in a separate variable. For example, 'BCleft=5'. In that way, they are easy to adapt and you prevent adding double fluxes.";//should become dict[1001]
//dict[2] = "You can use the nodes-vector to calculate the analytical solution for all x positions at the same time. Use substraction to look at the difference with the modelled heads. You can use the *dataframe.states()* function to get the modelled heads."; //should become dict[1002]
dict[1001] = "It can be handy to define the fluxes and boundary conditions in a separate variable. For example, 'BCleft=5'. In that way, they are easy to adapt and you prevent adding double fluxes.";//should become dict[1001]
dict[1002] = "You can use the nodes-vector to calculate the analytical solution for all x positions at the same time. Use substraction to look at the difference with the modelled heads. You can use the *dataframe.states()* function to get the modelled heads."; //should become dict[1002]

//dict[212] = "The model resembles the model made in Assignment 1 Part 2, except that now the external flux is now linearly decreasing with height.";
dict[2001] = "The model resembles the model made in Assignment 1 Part 2, except that now the external flux is now linearly decreasing with height.";
//dict[213] = "Check the documentation of the dataframe.* functions. The axes can be swapped by the argument named 'vertical', which takes a boolean.";
dict[2002] = "Check the documentation of the dataframe.* functions. The axes can be swapped by the argument named 'vertical', which takes a boolean.";
//dict[214] = "Check the documentation of dataframe.balance() to select an individual control volume.";
dict[2003] = "Check the documentation of dataframe.balance() to select an individual control volume.";
//dict[216] = "Read the documentation of the copy.model() function.";
dict[2004] = "Read the documentation of the copy.model() function.";
//dict[218] = "It might seem a bit odd to use the function oldstate(x) here since it is not defined yet. However, as long as we only define the storageflux function and not actually run it, we are allowed to use the undefined oldstate(x) function. We just have to make sure that the oldstate(x) function is defined and run before the storageflux calculation is run.";
dict[2005] = "It might seem a bit odd to use the function oldstate(x) here since it is not defined yet. However, as long as we only define the storageflux function and not actually run it, we are allowed to use the undefined oldstate(x) function. We just have to make sure that the oldstate(x) function is defined and run before the storageflux calculation is run.";
//dict[2120] = "Before the loop construct is entered, you can basically copy the previous code that describes the two time increments. Do not forget to assign the proper state values to the oldstate function before entering the while loop.";
dict[2006] = "Before the loop construct is entered, you can basically copy the previous code that describes the two time increments. Do not forget to assign the proper state values to the oldstate function before entering the while loop.";
//dict[2123] = "Create a function that simulates the potential temperature profile in time taking as its only argument the time step for which should be solved. Return a data container in which the profiles are stored. Run this function in a for loop for every time step and plot the result for comparison.";
dict[2007] = "Create a function that simulates the potential temperature profile in time taking as its only argument the time step for which should be solved. Return a data container in which the profiles are stored. Run this function in a for loop for every time step and plot the result for comparison.";
//dict[222] = "Think about the different sections in the model setup. First you define your system properties, then the mathematical part and finally the numerical part which makes the model ready to solve.";
dict[2008] = "Think about the different sections in the model setup. First you define your system properties, then the mathematical part and finally the numerical part which makes the model ready to solve.";
//dict[224] = "Several useful data frames are described in the documentation.";
dict[2009] = "Several useful data frames are described in the documentation.";
//dict[225] = "In order to add a graph to a figure, first the figure is created with the plot() function, thereafter lines() or points() can be used to add data.";
dict[2010] = "In order to add a graph to a figure, first the figure is created with the plot() function, thereafter lines() or points() can be used to add data.";
//dict[228] = "One should be able to initialize the states of the transient model by a combination of the functions do.initialize() and state.fun().";
dict[2011] = "One should be able to initialize the states of the transient model by a combination of the functions do.initialize() and state.fun().";
//dict[2211] = "What are the obligatory arguments in the system flux function? Which of these can be used for the saturated thickness?";
dict[2012] = "What are the obligatory arguments in the system flux function? Which of these can be used for the saturated thickness?";
//dict[2213] = "What is the time resolution of your recharge data?";
dict[2013] = "What is the time resolution of your recharge data?";
dict[3001] = "The concepts field capacity (pF 2) and wilting point (pF 4.2) are important here. When the water content of the soil is above the field capacity (pF<2) water will drain by gravity. The water in the soil below field capacity is the water that the soil can retain and supply to vegetation. However, when the soil water content is below wilting point (pF<4.2) the water is attached so strongly to the soil particles that plants cannot extract it any more. Therefore, the amount of water between field capacity and wilting point is considered to be the plant available water.";
dict[3002] ="To make the rainfun and evapfun available for your code chunk you can add 'source('topforcing.R')' as first line. You can also just run the script separately.";

//dict[41] = "The domain is defined by coordinates of two dimensions (x and y). Therefore, the domain is now defined by a 2 dimensional matrix. Type '?matrix()' in your console for more information.";
dict[4001] = "The domain is defined by coordinates of two dimensions (x and y). Therefore, the domain is now defined by a 2 dimensional matrix. Type '?matrix()' in your console for more information.";
//dict[42] = "Remember that a Neumann natural boundary condition is assigned to unspecified boundaries.";
dict[4002] = "Remember that a Neumann natural boundary condition is assigned to unspecified boundaries.";
//dict[45] = "Note that you can suppress graphing of the node generation by setting 'verbose=False'.";
dict[4003] = "Note that you can suppress graphing of the node generation by setting 'verbose=False'.";
//dict[413] = "Change both the domain values of 'xtopl' and 'ytopl' before interpolating the two dimensional model domain of TwoDmodel with the new transmissivities.";
dict[4004] = "Change both the domain values of 'xtopl' and 'ytopl' before interpolating the two dimensional model domain of TwoDmodel with the new transmissivities.";
//dict[417] = "add: welldist = func.dist.to.point(c(1000,1000),200,25,100) and set.discretisation(TwoDModelWell, nodes=list(type='geomgen', optdist=welldist, verbose=FALSE), method='FE') on separate lines to decrease the nodal distance around the location of the well.";
dict[4005] = "add: welldist = func.dist.to.point(c(1000,1000),200,25,100) and set.discretisation(TwoDModelWell, nodes=list(type='geomgen', optdist=welldist, verbose=FALSE), method='FE') on separate lines to decrease the nodal distance around the location of the well.";
//dict[419] = "Check the documentation of the func.dist.to.line() function in the FVFE2D package.";
dict[4006] = "Check the documentation of the func.dist.to.line() function in the FVFE2D package.";
//dict[420] = "Define domain range: plot(c(0,10), c(0,2), xlab='x', ylab='y', col='white') Plot polygon: polygon(domain, lwd=3, col=rgb(0.7,0.7,0.7,0.3)).";
dict[4007] = "Define domain range: plot(c(0,10), c(0,2), xlab='x', ylab='y', col='white') Plot polygon: polygon(domain, lwd=3, col=rgb(0.7,0.7,0.7,0.3)).";
//dict[422] = "Remember the function plot.BC.id().";
dict[4008] = "Remember the function plot.BC.id().";
//dict[425] = "The pressure function has 3 arguments: x, y and state. Its body contains one line: a return statement.";
dict[4009] = "The pressure function has 3 arguments: x, y and state. Its body contains one line: a return statement.";

// returns the hint for the corresponding question_number in the corresponding 
function showHint(question_number) {
	if (document.getElementById("Q"+question_number.toString()).innerHTML == " ") {
		var Qnum = question_number;
		document.getElementById("Q"+question_number.toString()).innerHTML = dict[Qnum];
	}
	else {
		var empty = " ";
		document.getElementById("Q"+question_number.toString()).innerHTML = empty;
	}
}

var coll = document.getElementsByClassName("collapsible");
var i;

for (i = 0; i < coll.length; i++) {
  coll[i].addEventListener("click", function() {
    this.classList.toggle("active");
    var content = this.nextElementSibling;
    if (content.style.display === "block") {
      content.style.display = "none";
    } else {
      content.style.display = "block";
    }
  });
} 