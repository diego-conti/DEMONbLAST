CURRENTDIR=$PWD

cd `dirname $0`
cd test
FLAGS="-Wctor-dtor-privacy -Wreorder -Wold-style-cast -Wsign-promo -Wchar-subscripts -Winit-self -Wmissing-braces -Wparentheses -Wreturn-type -Wswitch -Wtrigraphs -Wextra -Wno-sign-compare -Wno-narrowing -Wno-attributes -std=c++17"
EXTRA="../log.cpp"
for test in $(ls *.cpp); do
  executable=`basename $test .cpp`
  echo compiling $executable ...
  g++ -I/usr/local/include/wedge-0.3 $test $EXTRA -o $executable -O0 -g -lginac -lwedge -lcln -lgmp $FLAGS
  if [ $? -eq 0 ]; then ./$executable; fi  
done 

for out in $(ls *.out); do
  file_to_check=`basename $out .out`
  diff $file_to_check.out $file_to_check.ok
  if [ $? -eq 0 ]; then echo $file_to_check OK; 
  else echo $file_to_check ERROR; fi
done

cd $CURRENTDIR
