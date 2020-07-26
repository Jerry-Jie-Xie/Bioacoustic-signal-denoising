function MSE = calMSE(t1,t2)

error = t1-t2;               %Error between two signal
sqrd_err = error.^2;         %Taking element wise sum of the error vector
MSE = mean(sqrd_err,'all');  %use ¡®all¡¯ in mean function to mean over all the elements of
% squared error signal

end