# Lab Resolution

## Exercise 1
```math
I -> A, A, A (Data Parallelism)
A -> B, C (Functional Parallelism)
C -> A, A, A (Data Parallelism)
```

## Exercise 2
```math
6 = 1 / (f + (1 - f) / 10)
f = 2 / 27
TParallelizable = 1 - f = 25/27
```

## Exercise 3 (A)
```math
SUM[1 ≤ q ≤ p] [f(q, θ)] = 1
SUM[1 ≤ q ≤ p] [β * [(1 - θ) * (p + 1) + (2θ - 1) * q]] = 1
β * [p * (1 - θ) * (p + 1) + (2θ - 1) * SUM[1 ≤ q ≤ p] [q]] = 1
β * [p * (1 - θ) * (p + 1) + (2θ - 1) * p * (p + 1) / 2] = 1
β * p * (p + 1) * [1 - θ + (2θ - 1) / 2] = 1
β * p * (p + 1) * (1 - θ + θ - 1 / 2) = 1
β * p * (p + 1) * (1 / 2) = 1
β * p * (p + 1) / 2 = 1
β * p * (p + 1) = 2
β = 2 / [p * (p + 1)]
```

## Exercise 3 (B)
```math
S = 1 / SUM[1 ≤ q ≤ p] [f(q, θ) / q]
S = 1 / SUM[1 ≤ q ≤ p] [β * [(1 - θ) * (p + 1) + (2θ - 1) * q] / q]
S = 1 / [2 / [p * (p + 1)] * SUM[1 ≤ q ≤ p] [[(1 - θ) * (p + 1) + (2θ - 1) * q] / q]]
S = p * (p + 1) / [2 * SUM[1 ≤ q ≤ p] [[(1 - θ) * (p + 1) + (2θ - 1) * q] / q]]
S = p * (p + 1) / [2 * SUM[1 ≤ q ≤ p] [(1 - θ) * (p + 1) / q + (2θ - 1)]]
S = p * (p + 1) / [2 * (1 - θ) * (p + 1) * SUM[1 ≤ q ≤ p] [1 / q] + 2p * (2θ - 1)]
S = p * (p + 1) / [2 * (1 - θ) * (p + 1) * (ln(p) + 1) + 2p * (2θ - 1)]
```

## Exercise 3 (C)
```math
S = p * (p + 1) / [2 * (1 - θ) * (p + 1) * (ln(p) + 1) + 2p * (2θ - 1)]
For θ = 0:
S = p * (p + 1) / [2 * (p + 1) * (ln(p) + 1) - 2p]
For θ = 1 / 2:
S = p * (p + 1) / [(p + 1) * (ln(p) + 1)]
For θ = 1:
S = p * (p + 1) / 2p
```

## Exercise 4
```math
CPI = 2
OccupationTime = 1.25 * NI * (1 - 0.99)
TotalTime = CPI * NI = 2 * NI
OccupationRate = OccupationTime / TotalTime
OccupationRate = 1.25 * NI * (1 - 0.99) / [2 * NI]
OccupationRate = 1.25 * 0.01 / 2
OccupationRate = 0.00625
```

## Exercise 5 (A)
- Multiple Processors may access the same Pages, so the Common Cache can already have it loaded, without the need to access the Main Memory.
- There would be no replication of Pages through Multiple Individual Caches.
- Whenever a Write is performed there is no need to invalidate Multiple Individual Caches.

## Exercise 5 (B)
- Fewer Cache Swaps because each Process Cache only loads the addresses that the Process requires, so there may be less Cache Misses.
