# global parameters
################################################################################
VBase = 393.4#400
PBase = 1e4
IBase = PBase/VBase
resolution = 1e-3
samplerate = round(Int, 1/resolution)
################################################################################

# define structs for measurement data and complex phasor representation
################################################################################
mutable struct measurement
    t; U1; U2; U3; I1; I2; I3
end

mutable struct phasor
    t; u; uRot; uAmp; uPhi; i; iRot; iAmp; iPhi; P; Q; f; fAvg
end
################################################################################

# define transformations to complex phasor representation and post-processing
################################################################################
ABTransform = sqrt(2.0/3.0)*[1.0 -0.5 -0.5;
                             0.0 sqrt(3.0)*0.5 -sqrt(3.0)*0.5;
                             1.0/sqrt(2.0) 1.0/sqrt(2.0) 1.0/sqrt(2.0)]

function averageWin(data, winSize::Int)
    output = copy(data)
    for k in 1:length(data)-winSize
        output[k] = sum(data[k:k+winSize-1])/winSize
    end
    return output
end

function getFreqAvg(data, rate::Int)
    data_fft = rfft(data)
    val, fInd = findmax(abs2.(data_fft[2:end]))
    fSteps = rfftfreq(length(data), rate)
    fAvgOut = fSteps[fInd+1]
    return fAvgOut
end

function getDerivative(data, res)
    len = length(data)
    ddt1 = zeros(len)
    ddt3 = zeros(len)
    for k in 3:len-2
        ddt1[k] = (data[k+1] - data[k-1])/(2.0*res)
        ddt3[k] = (data[k+2]/2.0 - data[k+1] + data[k-1] - data[k-2]/2.0)/res^3
    end
    ddt1Out = (ddt1 .- res^2/6.0*ddt3)
    return ddt1Out
end

function getAngle(data)
    phiOut = angle.(data)
    for k in 2:length(data)
        count = 0
        while count < 1
            count += 1
            diff = abs(phiOut[k] - phiOut[k-1])
            if diff > 6.0
                count -= 1
                if phiOut[k] > phiOut[k-1]
                    phiOut[k] -= 2*pi
                else
                    phiOut[k] += 2*pi
                end
            end
        end
    end
    return phiOut
end

function toComplex(meas::measurement, phaseShift=0.)
    fAvg = getFreqAvg(Array{Float64}(meas.U1), samplerate)
    win = round(Int, 1/(resolution*fAvg))
    u_ab = ABTransform*transpose(hcat(meas.U1, meas.U2, meas.U3))
    u = (u_ab[1,:] .+ im*u_ab[2,:])/VBase.*exp.(im*phaseShift)
    uRot = u.*exp.(im*2*pi*fAvg*meas.t)
    uRotAvg = averageWin(uRot, win)
    uAmp = abs.(uRotAvg)
    uPhi = getAngle(uRotAvg)
    i_ab = ABTransform*transpose(hcat(meas.I1, meas.I2, meas.I3))
    i = (i_ab[1,:] .+ im*i_ab[2,:])/IBase
    iRot = i.*exp.(im*2*pi*fAvg*meas.t)
    iRotAvg = averageWin(iRot, win)
    iAmp = abs.(iRotAvg)
    iPhi = getAngle(iRotAvg)
    f = fAvg .- getDerivative(uPhi, resolution)/(2.0*pi)
    #S = uRotAvg.*conj.(iRotAvg)
    S = conj.(uRotAvg).*(iRotAvg)
    P = real.(S)
    Q = imag.(S)
    return phasor(meas.t, u, uRot, uAmp, uPhi, i, iRot, iAmp, iPhi, P, Q, f, fAvg)
end
################################################################################
