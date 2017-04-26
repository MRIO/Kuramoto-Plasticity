out1 = kuramotoSheet([10 10],1,'plasticity',[1 10 100])

close all

oscmean = mean(out1.oscillators);
oscstd = std(out1.oscillators);
find( (out1.oscillators > (oscmean+2.5*oscstd)) | (out1.oscillators < (oscmean-2.5*oscstd)))

imagesc(out1.connectivity)