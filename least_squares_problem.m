x=linspace(-1,1,200).';
y=x.^4-2*x.^3+x.^2-x+1;
matrixOfOnes = ones(200,1);
A=[matrixOfOnes x x.^2 x.^3]; %% do not use comma or colon here
leastSquaresCoeff = inv((A.'*A))*(A.')*y;
leastSquaresYValues = A*leastSquaresCoeff;
plot(x,y,'blue');
legend('','observed values');
xlabel('x values');
ylabel('y values');

hold
plot(x,leastSquaresYValues,'red'); %% least squares fit
legend('','leastSquares_y');

%% CORRECT SOLUTION