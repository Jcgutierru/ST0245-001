function draw_polyg_faces_2D(verts, loops, face_color)

low_pointer = 1;
N_faces = size(loops, 2);

for face=1:N_faces
    loop=face;
    up_pointer = low_pointer + loops(loop)-1;
    loop_i = verts(:,low_pointer:up_pointer);
    loop_i = [loop_i loop_i(:,1)];
    
    fill(loop_i(1,:),loop_i(2,:),face_color)
    
    hold on
    low_pointer = up_pointer+1;
end    





