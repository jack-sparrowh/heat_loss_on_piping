# heat_loss_on_piping
Code calculating heat loss (or gain) of liquid flowing at constant mass flow rate through insulated (or not) piping.

The code was prepared to run inside excels implementation of python IDE, that is meant to be utilized inside company I'm working for (who knows if it will). Because of that I had to use partials to load all the functions with the data after user specifies them.

The code is rough around the edges, because for some reason I convinced myself that writting it in functions only is a good idea. I'm rewritting it to utilize classes, it should be faster and more clean as I'm able to get rid of heaviside_sum function that does a lot of unnecessary calculations.

Most importantly this is a project that allowed me to understand steady-state heat transfer.

I'm adding also draft of my blog post where I'm going through all derivations. I'm writting all of it in my own free time, so it takes me a lot of time.

My only hope is that all the equations are rigorous under this list of assumptions:
	\begin{enumerate}
       	\item Steady-state operation (both mass and energy).
		\item The change in kinetic and potential energy is negligible.
		\item The flow is turbulent and pipes are full of liquid.
		\item The fluid in question is incompressible.
		\item The specific heat capacity is constant.
		\item No axial heat transfer.
		\item Piping is horizontal and enclosed by air only at ambient pressure.
	\end{enumerate}
