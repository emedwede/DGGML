#include <vtkActor.h>
#include <vtkDataSetMapper.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>

#include <vtkGlyph3DMapper.h>
#include <vtkSphereSource.h>
#include <vtkTubeFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkExtractEdges.h>
#include <vtkShrinkFilter.h>

int main(int argc, char* argv[])
{
  // parse command line arguments
  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " Filename(.vtu) e.g. graph.vtu"
              << std::endl;
    return EXIT_FAILURE;
  }

  std::string filename = argv[1];

  // read all the data from the file
  vtkNew<vtkXMLUnstructuredGridReader> reader;
  reader->SetFileName(filename.c_str());
  reader->Update();

  vtkNew<vtkNamedColors> colors;

  // Create a mapper and actor
  vtkNew<vtkDataSetMapper> mapper;
  mapper->SetInputConnection(reader->GetOutputPort());
  mapper->ScalarVisibilityOff();
  
  //extract the edges 
  vtkNew<vtkExtractEdges> extractEdges;
  extractEdges->SetInputConnection(reader->GetOutputPort());

  // Tube the edges
  vtkNew<vtkTubeFilter> tubes;
  tubes->SetInputConnection(extractEdges->GetOutputPort());
  tubes->SetRadius(0.1);
  tubes->SetNumberOfSides(21);
  vtkNew<vtkPolyDataMapper> edgeMapper;
  edgeMapper->SetInputConnection(tubes->GetOutputPort());
  edgeMapper->SetScalarRange(0, 26);
  vtkNew<vtkActor> edgeActor;
  edgeActor->SetMapper(edgeMapper);
  edgeActor->GetProperty()->SetSpecular(0.6);
  edgeActor->GetProperty()->SetSpecularPower(30);
  ;
 
  // Glyph the points
  vtkNew<vtkSphereSource> sphere;
  sphere->SetPhiResolution(21);
  sphere->SetThetaResolution(21);
  sphere->SetRadius(0.15);
  vtkNew<vtkGlyph3DMapper> pointMapper;
  pointMapper->SetInputConnection(reader->GetOutputPort());
  pointMapper->SetSourceConnection(sphere->GetOutputPort());
  pointMapper->ScalingOff();
  pointMapper->ScalarVisibilityOff();
  vtkNew<vtkActor> pointActor;
  pointActor->SetMapper(pointMapper);
  pointActor->GetProperty()->SetDiffuseColor(
      colors->GetColor3d("Banana").GetData());
  pointActor->GetProperty()->SetSpecular(0.6);
  pointActor->GetProperty()->SetSpecularColor(1.0, 1.0, 1.0);
  pointActor->GetProperty()->SetSpecularPower(100);
  ;

  vtkNew<vtkActor> actor;
  actor->SetMapper(mapper);
  actor->GetProperty()->EdgeVisibilityOn();
  actor->GetProperty()->SetLineWidth(2.0);
  actor->GetProperty()->SetColor(colors->GetColor3d("MistyRose").GetData());
  
  // The geometry
  vtkNew<vtkShrinkFilter> geometryShrink;
  geometryShrink->SetInputConnection(reader->GetOutputPort());
  geometryShrink->SetShrinkFactor(.8);

  vtkNew<vtkDataSetMapper> geometryMapper;
  geometryMapper->SetInputConnection(geometryShrink->GetOutputPort());
  geometryMapper->SetScalarModeToUseCellData();
  geometryMapper->SetScalarRange(0, 11);


  vtkNew<vtkActor> geometryActor;
  geometryActor->SetMapper(geometryMapper);
  geometryActor->GetProperty()->SetLineWidth(15);
  geometryActor->GetProperty()->EdgeVisibilityOn();
  geometryActor->GetProperty()->SetEdgeColor(0, 0, 0);


  vtkNew<vtkProperty> backFace;
  backFace->SetColor(colors->GetColor3d("Tomato").GetData());
  actor->SetBackfaceProperty(backFace);

  // Create a renderer, render window, and interactor
  vtkNew<vtkRenderer> renderer;
  vtkNew<vtkRenderWindow> renderWindow;
  renderWindow->AddRenderer(renderer);
  vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
  renderWindowInteractor->SetRenderWindow(renderWindow);
  
    
  // Add the actor to the scene
  renderer->AddActor(actor);
  renderer->AddActor(pointActor);
  renderer->AddActor(edgeActor);
  renderer->AddActor(geometryActor);
  renderer->SetBackground(colors->GetColor3d("Wheat").GetData());

  // Render and interact
  renderWindow->SetSize(640, 480);

  renderWindow->Render();
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}
