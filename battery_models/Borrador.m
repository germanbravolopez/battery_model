clear all
A = [];
for i=1:1:5
    A(1,i) = i;
end
for i=1:1:4
    A(1,5+i) = 0;
end
for i=1:1:5
    A(1,(9+(i)):(10+2*(i))) = [1 2];
end