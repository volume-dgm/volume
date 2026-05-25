
# volume-dgm
Программа для моделирования задач механики сплошных сред разрывным методом Галеркина.
## Содержание

- [Сборка проекта](#сборка-проекта)
- [Запуск программы](#запуск-программы)
- [Файлы конфигурации](#файлы-конфигурации)
  - [MeshBuilder](#meshbuilder)
  - [LambdaParser](#lambdaparser)
  - [Mesh](#mesh)
    - [MediumParams](#mediumparams)
    - [Boundaries](#boundaries)
    - [Contacts](#contacts)
  - [Snapshot](#snapshot)
  - [Schedule](#schedule)
  - [Task](#task)
  - [Solver](#solver)

## Сборка проекта

Осуществляется при помощи cmake. Зависимости, которые ожидаются установленными на системе:
- какая-нибудь реализация MPI
- OpenMP
- zlib

Зависимости, подгружаемые через git submodules:
- Metis (с GKlib)
- googletest


Для их загрузки нужно либо копировать репозиторий вместе с ними:
```
git clone --recurse-submodules -b af https://github.com/volume-dgm/volume.git
```
или инициализировать их после стандартного git clone: 
```
git submodule update --init --recursive.
```

Для сборки проекта необходимо стандартно создать папку build и вызвать из нее cmake и make (можно использовать любую другую привычную систему сборки, например, Ninja):
```
mkdir build && cd build 
cmake -DCMAKE_BUILD_TYPE=Release ..
make 
```
Принимаются следующие переменные:
- -DDEBUG_SYMBOLS=1: добавляет флаг компиляции -g (позволяет искать ошибки через gdb, но увеличивает размер бинарных файлов)
- -DCOMPILE_TESTS=1: включает сборку тестов через googletest (на данный момент пока недоделано и ничего на самом деле не тестирует, так что включать не обязательно)

## Запуск программы
Для запуска нужно предоставить программе размерность задачи, порядок используемых для аппроксимации решения полиномов, путь к файлу конфигурации (или просто его название) и расчетную сетку в правильном формате. Путь к последней указывается в конфиге, все остальное - передается либо как аргументы командной строки, либо как запись в файле task.xml. В первом случае запуск выглядит так:
```
./task dimCount polynomialsOrder taskName
```
- dimCount: размерность задачи (2 или 3)
- polynomialsOrder: порядок полиномов (число от 1 до 5; если необходимо другое - следует в файле src/Task/Task/Task.cpp добавить соответсвующую запись в switch-case на строках 30-55 и перекомпилировать программу)
- taskName:  либо путь к файлу конфигурации (определяется по наличию разрешения .xml), либо просто название конфига (в таком случае конфиг будет искаться файл config/$taskName.xml). То есть следующие два варианта аналогичны:
```
./task 2 1 test_diff
./task 2 1 config/test_diff.xml
```

При указании параметров в файле task.xml исполняемые файлы запускаются без аргументов:
```
./task
```
А внутри файла task.xml ожидается запись вида
```xml
<?xml version="1.0" encoding="UTF-8" ?>
<Settings fileName="config/test_diff.xml" dimsCount="2" polynomialsOrder="1"/>
```
Здесь аналогично в fileName можно передать название файла без расширения.

Для запуска исполняемого файла meshbuilder все аналогично, только не требуется указывать порядок полиномов:
```
./meshbuilder 2 test_diff
```
или
```
./meshbuilder 2 config/test_diff.xml
```
Внутри task.xml polynomialsOrder можно оставить - проигнорируется.

В программе реализовано распараллеливание с помощью OpenMP (запуск на систем с распределенной памятью через MPI пока не рассматриваем). По умолчанию задействуются все возможные потоки (можно узнать командой nproc на UNIX-подобных системах). Ограничивать число потоков можно через соответсвующие переменные среды в терминале:
```
export OMP_THREAD_LIMIT=1; ./task
```
## Файлы конфигурации 
Здесь со временем надеюсь описать все возможные параметры, пока же приведу отличия по сравнению с тем что было раньше для тех кто знает.

Для некоторых полей будет уточняться тип данных, к которому в конечном итоге приведется запись:
 - Scalar: вещественное число (e.g. "13.37")
 - unsigned: неотрицательное натуральное число (e.g. "6")
 - bool: "true" или "false"
 - Vector: вектор размерности dimsCount (e.g. "14.03 21.21" или "1 0 0"). zero-vector - вектор нулей.
 - Vector3: обязательно трехмерный вектор (e.g. ось вращения даже в 2d или "0 0 1" или "0 0 -1")
 - GenericDataType: [смотри здесь](#lambdaparser) - обобщенные скаляры (GenericdataType|Scalar) и векторы (GenericDataType|Vector). 
В случаях когда параметр является опциональным или имеется какое-либо значение по умолчанию рядом будет запись типа optional, default="". (optional == не ошибка не указать если не нужно, default == указать скорее всего надо, но если не будет, то применится данное значение)

### MeshBuilder

### LambdaParser

Из крупного - появился парсер лямбда-выражений в конфиге. По умолчанию на данный момент он работает с параметрами среды и начальными условиями; для контактных и граничных условий реализация пока не приводится потому что скорость работы программы в таком случае падает до катастрофически низких значений (возможно, оптимизируем в будущем, либо поймем, что это нам там и не надо вовсе).

Он работает по умолчанию для полей типа GenericDataType - для GenericDataType|Scalar можно указывать как лямбда-выражение (определяется по наличию подстроки lambda в начале), так и численное значение как раньше; для GenericDataType|Vector придется написать лямбда-выражение. В отдельном блоке в конфиге можно присвоить глобальные выражения для последующего использования командой AssignGlobal, которая принимает параметры name и value:
```xml
<LambdaParser> 
  <AssignGlobal name="it" value="5.0"/>
  <AssignGlobal name="ib" value="-5.0"/>
  <AssignGlobal name="Gmax" value="1.87e+9"/>
  <AssignGlobal name="Gmin" value="lambda global . / (get global 'Gmax') 2 "/>
  <AssignGlobal name="iniSpeed" value="2,0"/>
</LambdaParser>
```
- name: строка с именем global
- value: может быть как просто числом (принимаются double с точкой или запятой в качестве разделителя с поддержкой экспоненциальной записи), так и лямбда-выражением (в таком случае числа внутри тела выражения записываются только с запятой в качестве разделителя)

Затем в любом выражении применить присвоенный global можно при помощи функции get global 'globalName'.

В парсере по умолчанию определены global'ы для векторов:
```
"lambda global . 'Vector2'", "lambda global x y . insert (insert table 'x' x) 'y' y"
"lambda global . 'Vector3'", "lambda global x y z . insert (insert (insert table 'x' x) 'y' y) 'z' z"
```
Их следует использовать в полях, которые хранят в себе векторные величины; например, для скорости будет использована запись вида:
```
velocity="lambda global x y . (get global 'Vector2') (get global 'iniSpeed') 0"
```
В выражениях, которые возвращают данные уже непосредственно для использования в программе, следует указывать в качестве переменных (до точки) пространственные координаты: это обязательно для векторов и опционально для скалярных величин (в примерах выше двумерный вектор обязательно "lambda global x y . value_1 value_2"; скаляр G, не зависящий от координат, может быть как: 
```
"lambda global . / (get global 'Gmax') 2 "
```
, так и 
```
"lambda global x y . / (get global 'Gmax') 2 "
```
Для величин, которые меняются от координаты, соответственно, переменные указывать обязательно:
```
G="lambda global x y . (+ (get global 'Gmin') (/ (* (- y (get global 'ib')) (- (get global 'Gmax') (get global 'Gmin'))) (- (get global 'it') (get global 'ib'))))"
```

Поддерживается загрузка global'ов из внешнего файла. Это можно указать как в заголовке блока LambdaParser, так и командой LoadGlobals:
```xml
<LambdaParser loadFromFile="externalLambdas.xml">
  <LoadGlobals filename="someOtherFile.xml"/>
</LambdaParser>
```
Во внешних файлах ожидается аналогичный блок LambdaParser, который парсится точно так же как и в оригинальном конфиге.

### Mesh

Здесь задается путь к обработанной через meshbuilder сетке (потом мб поменяю тоже на поиск нужных файлов в стандартном месте через taskName) и несколько других геометрических параметров, в которых я пока, кс ожалению, не так разобрался, чтобы внятно объяснить, поэтому привожу старое их описание:
```xml
<Mesh fileName="meshes/3d_double_box[<domain>]" unfoldIterationsCount="0" minGridHeight="2e-4" moveMassCenter="false" collisionWidth= "0.01">
```
  - fileName: путь к файлам сетки без разрешения. \<domain> заменяется на индекс домена при разбиении сетки для вычислений на нескольких процессорах через MPI
  - unfoldIterationsCount: unsigned, optional, default="0" - number of iterations to try and fix collapsed mesh 
  - minGridHeight: Scalar, optional, default="0" - minimal cell height used by unfolding algorithm
  - moveMassCenter: bool, optional, default="false" - перемещает область snapshot'a вместе с сеткой (не проверял как работает)

  Затем следуют три дочерних элемента с параметрами сред, граничными и контактными условиями.

#### MediumParams

Теперь всегда используется PerSubmesh. Внутри себя он должен иметь столько дочерних элементов, сколько разных подсеток было задано в расчетной сетке:
```xml
<MediumParams>
  <PerSubmesh fileName="meshes/meshName.params">
    <Submesh index="0" ... parameters .../>
  </PerSubmesh>
</MediumParams>
```
fileName - путь к бинарному файлу с индексами используемого материала для каждой ячейки среды.

Внутри Submesh задаются следующие параметры:
 - Обязательно: пара упругих констант. Для любых из них используется GenericDataType|Scalar. Принимаются такие варианты:
  - lambda="2" mju="1": константы Ламе
  - E="2" nu="1": модуль Юнга и коэффициент Пуассона
  - E="2" G="1": модуль Юнга и модуль сдвига
  - pSpeed="2" sSpeed="1": скорости продольных и поперечных волн.
 - internalContactType="1": unsigned, default="0" - выставляемое между ячейками данного материала контактное условие
 - rho="920": GenericDataType|Scalar, default="1.0" - плотность
 - Опционально - параметры модели пластичности:
  - k="2.2e+5" alpha="0.1": GenericDataType|Scalar - параметры для расчета момента начала пластического течения
  - brittle="false": bool - потом разберусь, как-то то ли включает то ли выключает расчет пласчитности взависимости от параметров allowPlasticity и allowContinuousDestruction похоже
  - maxPlasticDeform="0.012": GenericDataType|Scalar - максимальная допустимая пластическая деформация
  - powderShearMult="0.1": Scalar, default="0" - тоже пока хз, по сути должно быть одной вещью но в коде кажись вся физика с ним связанная закоменчена и стоит что-то непонятное
 - fixed="false": bool, default="false" - фиксирована ли сетка

#### Boundaries

Используются следующие типы граничных условий:
 - Fixed 
 - Free 
 - Absorb
 
На каждый тип граничных условий, заданных в сетке, должна присутствовать соответствующая секция вида:
```xml
<Free interactionType="1" dynamicContactInteractionType = "2" />
```
  - interactionType: unsigned - индекс граничного условия в сетке
  - dynamicContactInteractionType: unsigned, optional, default="-1" - тип контактного условия, которое будет использовано при появлении динамического контакта с этой границей (значение -1 по умолчанию обозначает отсутствие динамического контакта (получается, не настолько unsigned, как казалось...))

На свободной границе можно задать внешнюю силу, на фиксированной - скорость. Делается это при помощи векторных функторов:
```xml
<Fixed interactionType="0" dynamicContactInteractionType = "2" >
  <ExternalVelocity>
    <RotatingFunctor pos="-5 0" rotationAxis="0 0 1" angularVelocity="5" linearVelocity="2 0" linearTime="10"/>
  </ExternalVelocity>
</Fixed>
```
Для силы:
```xml
<Free interactionType= "0" dynamicContactInteractionType= "1"> 
  <ExternalForce> 
    <ConstFunctor value= "10 0 0"/> 
  </ExternalForce> 
<Free>
```
Таких функторов может быть и несколько для одной границы. Они бывают следующих типов (сейчас пока перечислены только те, с которыми лично я имел дело):

  - ConstFunctor - постоянное значение:
  - ```xml
    \<ConstFunctor value= "10 0 0"/>
    ```
    - value: Vector - значение скорости/силы
  - RotatingFunctor - вращение вокруг оси:
     ```xml
       <RotatingFunctor pos="-5 0" rotationAxis="0 0 1" angularVelocity="5" linearVelocity="2 0" linearTime="10"/>
      ```
    - pos: Vector - радиус-вектор точки на оси вращения
    - rotationAxis: Vector3 - направление оси вращения 
    - angularVelocity: Scalar - угловая скорость (рад/с)
    - linearVelocity: Vector, optional, default=zero-vector - линейная скорость движения границы (для комбинированного поступательного и вращательного движения)
    - linearTime: Scalar, optional, default="-1" - продолжительность поступательного движения (если currTime > linearTime то считается что linearVelocity=0; если linearTime=-1 то linearVelocity работает на протяжении всего времени моделирования)

Тут стоит учитывать, что внешняя скорость задается только на граничных элементах, а не по всему телу. Поэтому при начале движения, если не задать те же скорости как начальное условие по области, будут возникать начальные напряжения в теле вплоть до разрушения при больших значениях скорости. Как для вращения, так и для поступательного движения есть соответствующие начальные условия - используйте их.

#### Contacts

### Snapshot

Здесь перечисляется когда и какие файлы сохранять. Весь вывод теперь пишется в директорию out/taskName, если ее нет - она создается, если она есть - ее содержимое на всякий случай копируется в out/taskName_backup. Все сохраянемые файлы имеют имя вида "\<fileType>[<domain>]_<step>.vtk (.vti for data)", где
  - \<fileType> - тип записанных данных (т.е или data, или mesh, или cellInfos/contacts/boundaries)
  - \<domain> - номер вычислительного узла (всегда будет 000 без использования MPI)
  - <step> - номер итерации

Таким образом, имена файлов больше писать не нужно и весь блок имеет вид такой:
```xml
<Snapshot>
  <Period time="5e-4" />
  <Data  writeVelocity="True" writeTension="True">
    <Box boxPoint1="-10.1 -5.0" boxPoint2="10.1 5.0" />
    <Resolution resolution="402 200" />
  </Data>
  <Mesh/>
  <CellInfos/>
</Snapshot>
```
В Period можно использовать или time="5e-4": Scalar, или frames="100": unsigned для сохранение каждые 5e-4 с или 100 итераций соответственно. 

остальное потом

### Schedule

### Task

Здесь указывается конечное время моделирования и начальные условия:
```xml
<Task destinationTime = "0.05">
  <IniState>
    <BoxState ... />
  </IniState>
</Task>
```

IniState бывают следующих видов (опять же, пока что только те, с которыми я имел дело, в дальнейшем протестирую все и добавлю):
- BoxState - начальная скорость в прямоугольной области:
```xml
<BoxState boxPoint1="0.0 -5.0" boxPoint2="10.1 5.0" submeshToApply="1" velocity="lambda global x y . 0 0"/>
```
  - boxPoint1, boxPoint2: Vector - крайние точки области, внутри которой будет применяться начальное условие
  - submeshToApply: unsigned, default=0: индекс подсетки (т.е значение из Submesh index="" в секции MediumParams), к которому будет применяться данное условие
  - velocity: GenericDataType|Vector: векторное лямбда-выражение со значением скорости
- BoxRotationState - начальная скорость при вращении вокруг оси:
```xml
<BoxRotationState boxPoint1 = "-100.0 -50.0" boxPoint2 = "100.1 50.0" submeshToApply="0" pos="-5 0" rotationAxis="0 0 1" angularVelocity="5" linearVelocity="2 0" />
```
  - boxPoint1, boxPoint2: Vector - крайние точки области, внутри которой будет применяться начальное условие
  - submeshToApply: unsigned, default=0: индекс подсетки (т.е значение из Submesh index="" в секции MediumParams), к которому будет применяться данное условие
  - pos: Vector - радиус-вектор точки на оси вращения
  - rotationAxis: Vector3 - направление оси вращения 
  - angularVelocity: Scalar - угловая скорость (рад/с)
  - linearVelocity: Vector - линейная скорость движения (для комбинированного поступательного и вращательного движения)

### Solver
