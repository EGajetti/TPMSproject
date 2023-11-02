#include "Utils.h"


// Implementation of the function defined in Utils.h

void printTime(clock_t start, clock_t end) {
	double elapsed = double(end - start) / CLOCKS_PER_SEC;
	cout << "Time measured: " << setprecision(7) << elapsed << " seconds.\n";
}


string* readConfiguration(const char* filename) {
	ifstream file(filename);

	string line;
	string var_name[15];
	string* var_value = new string[15];
	string inSection;
	size_t posequal;
	int var_pos = 0;

	while (getline(file, line)) {
		if (!line.length()) continue;
		if (line[0] == '#') continue;
		posequal = line.find('=');
		var_name[var_pos] = line.substr(0, posequal);
		if (posequal == line.size())
			var_value[var_pos] = "";
		else
			var_value[var_pos] = line.substr(posequal + 1, line.size());
		var_value[var_pos].erase(remove_if(var_value[var_pos].begin(), var_value[var_pos].end(), ::isspace), var_value[var_pos].end());

		var_pos++;
	}

	for (int i = 0; i < var_pos; i++)
		cout << var_name[i] << " ----> " << var_value[i] << endl;
	return var_value;
}


double* convertOrigin(string origin) {
	int i = 1;
	string x, y, z;
	static double orig[3];
	while (origin[i] != ',') {
		x.push_back(origin[i]);
		i++;
	}
	orig[0] = stod(x);
	i++;
	while (origin[i] != ',') {
		y.push_back(origin[i]);
		i++;
	}
	orig[1] = stod(y);
	i++;
	while (origin[i] != ']') {
		z.push_back(origin[i]);
		i++;
	}
	orig[2] = stod(z);
	return orig;
}


#ifdef GRAPHICAL

void renderSurface(vtkFlyingEdges3D* surface, vtkQuadricDecimation* decimate) {

	vtkNew<vtkNamedColors> colors;

	vtkNew<vtkRenderer> renderer;
	renderer->SetBackground(colors->GetColor3d("DarkSlateGray").GetData());

	vtkNew<vtkRenderWindow> renderWindow;
	renderWindow->AddRenderer(renderer);
	renderWindow->SetWindowName("MarchingCubes");

	vtkNew<vtkRenderWindowInteractor> interactor;
	interactor->SetRenderWindow(renderWindow);

	vtkNew<vtkPolyDataMapper> mapper;
	mapper->SetInputConnection(decimate->GetOutputPort());

	mapper->ScalarVisibilityOff();

	vtkNew<vtkActor> actor;
	actor->SetMapper(mapper);
	actor->GetProperty()->SetColor(colors->GetColor3d("MistyRose").GetData());

	renderer->AddActor(actor);

	renderWindow->Render();
	interactor->Start();
}
#endif // GRAPHICAL

