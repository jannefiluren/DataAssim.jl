
# Time correlated noise

function timecorr_noise(q, timestep, timedecorr)

	alfa = 1 - timestep / timedecorr

	q = alfa*q + sqrt(1-alfa^2)*randn()

	return q

end


# Resampling function

function resample(wk)

	Ns = length(wk)

	u = cumprod(rand(Ns) .^ (1./collect(Float64, Ns:-1:1)))

	u = u[length(u):-1:1]

	wc = cumsum(wk)

	label = zeros(Int64, Ns)

	k = 1
	for i = 1:Ns
		while wc[k] < u[i]
			k = k + 1
		end
		label[i] = k
	end

	return label

end

