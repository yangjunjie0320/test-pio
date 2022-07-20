import os
import sys

VMD_PATH = "/Users/yangjunjie/Applications/VMD\ 1.9.4a48-Catalina-Rev7.app/Contents/vmd/vmd_MACOSXX86_64"

def show_cube(inp_file, out_file="test.png", iso="0.08"):
    vmd_cmd = '''
    proc show_cube {filename {isoval 0.05}} {
    #default material
    set mater Glossy

    color Display Background white
    display depthcue off
    display rendermode GLSL
    axes location Off
    color Name C tan
    color change rgb tan 0.700000 0.560000 0.360000
    material change mirror Opaque 0.15
    material change outline Opaque 4.000000
    material change outlinewidth Opaque 0.5
    material change ambient Glossy 0.1
    material change diffuse Glossy 0.600000
    material change opacity Glossy 0.75
    material change shininess Glossy 1.0
    light 3 on

    foreach i [molinfo list] {
    mol delete $i
    }

    mol new $filename
    mol modstyle 0 top CPK 0.800000 0.300000 22.000000 22.000000
    mol addrep top
    mol modstyle 1 top Isosurface $isoval 0 0 0 1 1
    mol modcolor 1 top ColorID 12
    mol modmaterial 1 top $mater
    mol addrep top
    mol modstyle 2 top Isosurface -$isoval 0 0 0 1 1
    mol modcolor 2 top ColorID 22
    mol modmaterial 2 top $mater
    display distance -2.0
    display height 10
    }

    proc cubiso {isoval} {
    mol modstyle 1 top Isosurface $isoval 0 0 0 1 1
    mol modstyle 2 top Isosurface -$isoval 0 0 0 1 1
    }

    proc cub2 {filename1 filename2 {isoval 0.05}} {
    #default material
    set mater Glossy

    color Display Background white
    display depthcue off
    display rendermode GLSL
    axes location Off
    color Name C tan
    color change rgb tan 0.700000 0.560000 0.360000
    material change mirror Opaque 0.15
    material change outline Opaque 4.000000
    material change outlinewidth Opaque 0.5
    material change ambient Glossy 0.1
    material change diffuse Glossy 0.600000
    material change opacity Glossy 0.75
    material change shininess Glossy 1.0

    light 0 on
    light 1 on
    light 2 on
    light 3 on

    foreach i [molinfo list] {
    mol delete $i
    }

    mol new $filename1
    mol modstyle 0 top CPK 0.800000 0.300000 22.000000 22.000000
    mol addrep top
    mol modstyle 1 top Isosurface $isoval 0 0 0 1 1
    mol modcolor 1 top ColorID 12
    mol modmaterial 1 top $mater

    mol new $filename2
    mol modstyle 0 top CPK 0.800000 0.300000 22.000000 22.000000
    mol addrep top
    mol modstyle 1 top Isosurface $isoval 0 0 0 1 1
    mol modcolor 1 top ColorID 22
    mol modmaterial 1 top $mater

    display distance -2.0
    display height 10
    }

    proc cub2iso {isoval} {
    foreach i [molinfo list] {
    mol modstyle 1 $i Isosurface $isoval 0 0 0 1 1
    }
    }

    show_cube {%s} %6.4f
    rotate x to -70 
    render TachyonInternal {%s}
    '''%(inp_file, iso, out_file)

    with open("tmp.vmd", "w") as f:
        f.write(vmd_cmd)

    os.system(f"{VMD_PATH} < tmp.vmd")

if __name__ == "__main__":
    inp_file = sys.argv[1]
    out_file = sys.argv[2]
    iso      = float(sys.argv[3])

    show_cube(inp_file, out_file=out_file, iso=iso)
