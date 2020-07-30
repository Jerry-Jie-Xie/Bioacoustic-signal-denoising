% calculate signal to noise ratio
function SNR = GetSNR(signal, errorSignal)


type = 'no_sqrt';
switch type   
    
    case 'sqrt'
        SNR = 10 * log10(sqrt(mean(signal.^2))/sqrt(mean(errorSignal.^2)));
        
    case 'no_sqrt'
        SNR = 10 * log10(mean(signal.^2)/mean(errorSignal.^2));        
end

end
%[EOF]