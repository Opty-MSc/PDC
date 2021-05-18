RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'

# Navigate to Program Version
cd "../MPI" || exit

# Compile Program
make

# Run Test 'ex-2-5-0'
echo -e "${NC}ex-2-5-0"
OMP_NUM_THREADS=${1} mpirun -np ${2} ballAlg.o 2 5 0 >ex-2-5-0.tree
./ballQuery.o ex-2-5-0.tree 3 1 >ex-2-5-0.query

if cmp --silent -- ex-2-5-0.query ../Tests/Expected/ex-2-5-0.query; then
  echo -e "${GREEN}Success!"
else
  echo -e "${RED}Failed!"
fi

# Run Test 'ex-2-8-0'
echo -e "${NC}ex-2-8-0"
OMP_NUM_THREADS=${1} mpirun -np ${2} ballAlg.o 2 8 0 >ex-2-8-0.tree
./ballQuery.o ex-2-8-0.tree 8 8 >ex-2-8-0.query

if cmp --silent -- ex-2-8-0.query ../Tests/Expected/ex-2-8-0.query; then
  echo -e "${GREEN}Success!"
else
  echo -e "${RED}Failed!"
fi

# Run Test 'ex-20-1000000-0'
echo -e "${NC}ex-20-1000000-0"
OMP_NUM_THREADS=${1} mpirun -np ${2} ballAlg.o 20 1000000 0 >ex-20-1000000-0.tree
./ballQuery.o ex-20-1000000-0.tree 1 2 3 4 5 6 7 8 9 1 2 3 4 5 6 7 8 9 1 2 >ex-20-1000000-0.query

if cmp --silent -- ex-20-1000000-0.query ../Tests/Expected/ex-20-1000000-0.query; then
  echo -e "${GREEN}Success!"
else
  echo -e "${RED}Failed!"
fi

# Run Test 'ex-3-5000000-0'
echo -e "${NC}ex-3-5000000-0"
OMP_NUM_THREADS=${1} mpirun -np ${2} ballAlg.o 3 5000000 0 >ex-3-5000000-0.tree
./ballQuery.o ex-3-5000000-0.tree 4 5 6 >ex-3-5000000-0.query

if cmp --silent -- ex-3-5000000-0.query ../Tests/Expected/ex-3-5000000-0.query; then
  echo -e "${GREEN}Success!"
else
  echo -e "${RED}Failed!"
fi

# Run Test 'ex-4-10000000-0'
echo -e "${NC}ex-4-10000000-0"
OMP_NUM_THREADS=${1} mpirun -np ${2} ballAlg.o 4 10000000 0 >ex-4-10000000-0.tree
./ballQuery.o ex-4-10000000-0.tree 2 4 6 8 >ex-4-10000000-0.query

if cmp --silent -- ex-4-10000000-0.query ../Tests/Expected/ex-4-10000000-0.query; then
  echo -e "${GREEN}Success!"
else
  echo -e "${RED}Failed!"
fi

# Run Test 'ex-3-20000000-0'
echo -e "${NC}ex-3-20000000-0"
OMP_NUM_THREADS=${1} mpirun -np ${2} ballAlg.o 3 20000000 0 >ex-3-20000000-0.tree
./ballQuery.o ex-3-20000000-0.tree 1 5 9 >ex-3-20000000-0.query

if cmp --silent -- ex-3-20000000-0.query ../Tests/Expected/ex-3-20000000-0.query; then
  echo -e "${GREEN}Success!"
else
  echo -e "${RED}Failed!"
fi

# Run Test 'ex-4-20000000-0'
echo -e "${NC}ex-4-20000000-0"
OMP_NUM_THREADS=${1} mpirun -np ${2} ballAlg.o 4 20000000 0 >ex-4-20000000-0.tree
./ballQuery.o ex-4-20000000-0.tree 8 6 4 2 >ex-4-20000000-0.query

if cmp --silent -- ex-4-20000000-0.query ../Tests/Expected/ex-4-20000000-0.query; then
  echo -e "${GREEN}Success!"
else
  echo -e "${RED}Failed!"
fi

# Cleanup Program
rm -f ./*.tree ./*.query
make clean

# Navigate Back
cd ../Tests || exit
