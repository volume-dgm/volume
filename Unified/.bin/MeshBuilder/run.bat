if EXIST "meshes" (rmdir /wait /s /q meshes)
md "meshes"
start /wait MeshBuilder.exe > log.txt
xcopy ".\meshes\*.params" ".\..\Task\meshes" /s /y /q
xcopy ".\meshes\*.mesh" ".\..\Task\meshes" /s /y /q
