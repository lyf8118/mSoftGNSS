function llrValues = initLlr(receivedSignal, noiseVariance)
%INITLLR Calculate BPSK LLR values for AWGN soft-decision decoding.

noiseVariance = max(noiseVariance, eps);
channelReliability = 2 / noiseVariance;
llrValues = channelReliability * receivedSignal;
end
