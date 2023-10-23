function dvardt = varEvolution(t,var,s)

[time,varianceDissipation] = microMixing(s);
varDiss=interp1(time,varianceDissipation,t);
dvardt=varDiss*var;
end