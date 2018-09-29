
clc;
close all;
clear all;

fidO = fopen('output_costs.txt','w');
fidE = fopen('output_numiters.txt','w');


for Gamma = 0:1

fid = fopen('input_1.txt');
TotalVertices = fgetl(fid);
TotalVertices = str2num(TotalVertices);
StartingVertex = fgetl(fid);
StartingVertex = str2num(StartingVertex);
GoalVertex = fgetl(fid);
GoalVertex = str2num(GoalVertex);

V(1:TotalVertices,1) = Inf;%high-default value
V(StartingVertex) = 0;

B = zeros(TotalVertices, 1);

MainArray = dlmread('input_1.txt',' ',3,0);
Coordinates = dlmread('coords_1.txt',' ',0,0);
fclose(fid);

[V, B, Closedlist, Openlist] = DijkstraAStar(V,TotalVertices,StartingVertex,GoalVertex,Gamma, B, MainArray, Coordinates);

fidO = fopen('output_costs.txt','a');
fprintf(fidO,'%f ', V(GoalVertex));
fclose(fidO);

fidE = fopen('output_numiters.txt','a');
sizeOfClosedList = size(Closedlist);
fprintf(fidE,'%d ', sizeOfClosedList(1));
fclose(fidE);
end

fidE = fopen('output_numiters.txt','a');
fprintf(fidE,'\n');
fclose(fidE);

fidO = fopen('output_costs.txt','a');
fprintf(fidO,'\n');
fclose(fidO);

for Gamma = 0:1

fid = fopen('input_2.txt');
TotalVertices = fgetl(fid);
TotalVertices = str2num(TotalVertices);
StartingVertex = fgetl(fid);
StartingVertex = str2num(StartingVertex);
GoalVertex = fgetl(fid);
GoalVertex = str2num(GoalVertex);

V(1:TotalVertices,1) = Inf;%high-default value
V(StartingVertex) = 0;

B = zeros(TotalVertices, 1);

MainArray = dlmread('input_2.txt',' ',3,0);
Coordinates = dlmread('coords_2.txt',' ',0,0);
fclose(fid);

[V, B, Closedlist, Openlist] = DijkstraAStar(V,TotalVertices,StartingVertex,GoalVertex,Gamma, B, MainArray, Coordinates);


fidO = fopen('output_costs.txt','a');
fprintf(fidO,'%f ', V(GoalVertex));
fclose(fidO);

fidE = fopen('output_numiters.txt','a');
sizeOfClosedList = size(Closedlist);
fprintf(fidE,'%d ', sizeOfClosedList(1));
fclose(fidE);
end

fidE = fopen('output_numiters.txt','a');
fprintf(fidE,'\n');
fclose(fidE);

fidO = fopen('output_costs.txt','a');
fprintf(fidO,'\n');
fclose(fidO);

for Gamma = 0:1

fid = fopen('input_3.txt');
TotalVertices = fgetl(fid);
TotalVertices = str2num(TotalVertices);
StartingVertex = fgetl(fid);
StartingVertex = str2num(StartingVertex);
GoalVertex = fgetl(fid);
GoalVertex = str2num(GoalVertex);

V(1:TotalVertices,1) = Inf;%high-default value
V(StartingVertex) = 0;

B = zeros(TotalVertices, 1);

MainArray = dlmread('input_3.txt',' ',3,0);
Coordinates = dlmread('coords_3.txt',' ',0,0);
fclose(fid);

[V, B, Closedlist, Openlist] = DijkstraAStar(V,TotalVertices,StartingVertex,GoalVertex,Gamma, B, MainArray, Coordinates);

fidO = fopen('output_costs.txt','a');
fprintf(fidO,'%f ', V(GoalVertex));
fclose(fidO);

fidE = fopen('output_numiters.txt','a');
sizeOfClosedList = size(Closedlist);
fprintf(fidE,'%d ', sizeOfClosedList(1));
fclose(fidE);
end

function[V, B, Closedlist, Openlist] = DijkstraAStar(V,TotalVertices,StartingVertex,GoalVertex,Gamma, B, MainArray, Coordinates)

Openlist= {StartingVertex};
Closedlist = {};

V(1:TotalVertices,1) = Inf;
V(StartingVertex,1) = 0;

currentvertex = StartingVertex;

Finished = false;
while ~Finished
    
    CostToCome={};
    for pos = 1:size(Openlist)
        tempo = size(CostToCome);
        
        OpenVertex = Openlist{pos};
        
        if isempty(OpenVertex)
            continue;
        end
        
        %         if  Openlist{pos} == currentvertex
        %             CostToCome{tempo(1)+1,1} = 0;
        %         else
        %             TempRowNum = find(MainArray(:,1) == OpenVertex & MainArray(:,2) == currentvertex);
        %             CostToCome{tempo(1)+1,1} = MainArray(TempRowNum,3) + V(Openlist{pos});
        CostToCome{tempo(1)+1,1} = V(Openlist{pos}) + (Gamma * H(Coordinates,GoalVertex,Openlist{pos}));
        %         end
        
    end
    
    [MinVal, MinIndex] = min(cell2mat(CostToCome));
    CTCSize = size(CostToCome);
    if MinIndex < 1 || MinIndex > CTCSize(1)
        MinIndex = 1;
    end
    
    closedsize = size(Closedlist);
    Closedlist{closedsize(1)+1,1} = Openlist{MinIndex,1};
    currentvertex = Openlist{MinIndex,1};
    Openlist(MinIndex,:) = [];
    
    if currentvertex ~=GoalVertex
        NeighbourRowIndices = find(MainArray(:,2) == currentvertex);
        NumNeighbours = size(NeighbourRowIndices);
        for index = 1:NumNeighbours
            Neighbour = MainArray(NeighbourRowIndices(index),1);
            
            closedsize = size(Closedlist);
            PresentInClosedList = false;
            if find([Closedlist{:}] == Neighbour)
                PresentInClosedList = true;
                continue;
            end
            
            opensize = size(Openlist);
            PresentInOpenList = false;
            if find([Openlist{:}] == Neighbour)
                PresentInOpenList = true;
                continue;
            end
            
            Size = size(Openlist);
            Openlist{Size(1)+1,1}= Neighbour;
            
            Vnew = MainArray(NeighbourRowIndices(index),3) + V(currentvertex);
            
            if V(Neighbour) > Vnew
                V(Neighbour) = Vnew;
            end
            B(Neighbour) = currentvertex;
            
        end
        
    else
        Finished = true;
    end
    
end
end

function[distance] = H(Coordinates, GoalVertex, GivenVertex)
x1 = Coordinates(GoalVertex,1);
y1 = Coordinates(GoalVertex,2);
x2 = Coordinates(GivenVertex,1);
y2 = Coordinates(GivenVertex,2);

distance = sqrt((x2-x1)^2+(y2-y1)^2);
end


