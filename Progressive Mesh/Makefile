
compile = g++ -Wno-deprecated -pthread

all:	run

clean:
	rm ProgressiveMesh

run:	main.cpp
	$(compile) main.cpp Render.cpp PM.cpp MeshLib_Core/Mesh.cpp MeshLib_Core/Vertex.cpp -o ProgressiveMesh -lGL -lGLU -lglut 
