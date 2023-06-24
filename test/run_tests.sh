
all_tests=$(ls test-*.py)

echo "Running all the test-*.py to test the implementation"
echo ${all_tests}
echo "-----------"


for t in ${all_tests}; do
  echo ${t}
  python3 ${t}
done

