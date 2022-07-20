for f in $(cat cube_list.log); do
    echo "Processing $f"
    python show_cube.py ./cube/$f.cube ./png/$f.png 0.10
done
