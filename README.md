
# See documentation at [this site](http://brianinvo.github.io/xrayr/index.html)

# see [this](https://medium.com/@sebagomez/installing-the-docker-client-on-ubuntus-windows-subsystem-for-linux-612b392a44c4
# ) 
# for windows docker install instrux .... needed for neuralenhance



# neuralenahance

docker run --rm -v `pwd`:/ne/input -it alexjc/neural-enhance --help

alias enhance='function ne() { docker run --rm -v "$(pwd)/`dirname ${@:$#}`":/ne/input -it alexjc/neural-enhance ${@:1:$#-1} "input/`basename ${@:$#}`"; }; ne'

# Now run any of the examples above using this alias, without the `.py` extension.
# enhance --zoom=1 --model=repair images/broken.jpg

