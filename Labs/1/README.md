# Lab Resolution

## Exercise 1
```math
FractionOfTime in C = 0.0002
FractionOfTime in B = R
FractionOfTime in A = 1 - 0.0002 - R

S = 80
1 / (0.0002 + R / 50 + (1 - 0.0002 - R) / 100) = 80
0.0002 + R / 50 + (1 - 0.0002 - R) / 100 = 1 / 80
R / 50 + (1 - 0.0002 - R) / 100 = 1 / 80 - 0.0002
R / 50 + 1 / 100 - 0.0002 / 100 - R / 100 = 1 / 80 - 0.0002
2 * R / 100 - R / 100 = 1 / 80 - 0.0002 - 1 / 100 - 0.0002 / 100
R / 100 = 1 / 80 - 0.0002 - 1 / 100 - 0.0002 / 100
R = 100 * (1 / 80 - 0.0002 - 1 / 100 - 0.0002 / 100)
R ~= 0.23
```

## Exercise 2
```math
1GHz = 10^9Hz
ClockCycleTime = 1 / 10^9
RemoteMemoryAccessTime = 400 * 10^-9 Seconds

I = ApplicationInstructions
TimeWithoutRemoteMemoryAccess = I / 2 * ClockCycleTime
TimeWithRemoteMemoryAccess = 0.002 * I * RemoteMemoryAccessTime + (1 - 0.002) * I / 2 * ClockCycleTime

S = TimeWithRemoteMemoryAccess / TimeWithoutRemoteMemoryAccess
S = (0.002 * I * 400 * 10^-9 + (1 - 0.002) * I / 2 * 1 / 10^9) / (I / 2 * 1 / 10^9)
S = (0.002 * 1 * 400 * 10^-9 + (1 - 0.002) * 1 / 2 * 1 / 10^9) / (1 / 2 * 1 / 10^9)
S ~= 2.60
```

## Exercise 3 (A)
```math
[StartUp/Finish]Time = 10^-6 Seconds
IntTransferTime = 4 * 10^-9 Seconds
IterationTime = 100 * 10^-9 Seconds
MultiProcessTime = 2 * [StartUp/Finish]Time + 2 * IntTransferTime + SIZE / N * (IterationTime + IntTransferTime)
SingleProcessTime = 2 * [StartUp/Finish]Time + SIZE * IterationTime

S(N, SIZE) = SingleProcessTime / MultiProcessTime
S(N, SIZE) = (2 * 10^-6 + SIZE * 100 * 10^-9) / (2 * 10^-6 + 2 * 4 * 10^-9 + SIZE / N * (100 * 10^-9 + 4 * 10^-9))
```

## Exercise 3 (B)
```math
S(100, 100) = (2 * 10^-6 + 100 * 100 * 10^-9) / (2 * 10^-6 + 2 * 4 * 10^-9 + 100 / 100 * (100 * 10^-9 + 4 * 10^-9))
S(100, 100) ~= 5.68
S(100, 10000) = (2 * 10^-6 + 10000 * 100 * 10^-9) / (2 * 10^-6 + 2 * 4 * 10^-9 + 10000 / 100 * (100 * 10^-9 + 4 * 10^-9))
S(100, 10000) ~= 80.75
```
