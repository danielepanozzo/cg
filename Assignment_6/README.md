Rasterization
=============

In this exercise you will implement an interactive 2D editor for vector graphics. The editor will allow to draw simple shapes interactively.

### Using Eigen

In all exercises you will need to do operations with vectors and matrices. To simplify the code, you will use [Eigen](http://eigen.tuxfamily.org/).
Have a look at the [Getting Started](http://eigen.tuxfamily.org/dox/GettingStarted.html) page of Eigen as well as the [Quick Reference](http://eigen.tuxfamily.org/dox/group__QuickRefPage.html}) page for a reference of the basic matrix operations supported.

### SDL2

To handle input events (keyboard and mouse), we rely on the SDL2 library to abstract the low-level OS-specific details of the process. Please refer to the [documentation](https://www.libsdl.org) of SDL2 for more details.

### Preparing the Environment and Submission

Follow the instructions on the [general instructions page](../RULES.md) to set up what you need for the assignment.

Ex.1: Triangle Soup Editor
--------------------------

Implement an interactive application that allows to add, edit, and delete triangles. The following operations should be supported:

- The key <kbd>i</kbd> will enable triangle insertion mode. When this mode is enabled, every triple of subsequent mouse clicks will create a triangle connecting the three locations where the mouse has been pressed. The first click will create the starting point of the segment, which will be immediately visualized. As the mouse is moved, a preview of a segment will appear. After the second mouse click, a preview of the triangle will appear. After the third click, the current preview will transform into the final triangle.

![image](img/i.png)

- The key <kbd>o</kbd> will enable triangle translation mode. Each mouse click will select the triangle below the cursor (which will be highlighted), and every movement of the mouse (while keeping the button pressed) will result in a corresponding translation of the triangle. Note that the triangle should move on screen by the same amount as the cursor.

![image](img/o.png)

- The key <kbd>p</kbd> will enable delete mode. When the mouse is pressed, the triangle below the cursor is deleted.


Ex.2: Rotation/Scale
--------------------

When triangle translation mode is enabled, keep the current primitive selected after the mouse is released. If a primitive is selected and you press the keys <kbd>h</kbd> and <kbd>j</kbd>, the triangle will rotate by 10 degree clockwise or counter-clockwise, respectively. The rotations should be done around its barycenter, i.e. the barycenter of the triangle should not change. When the keys <kbd>k</kbd> or <kbd>l</kbd> are pressed, the primitive should be scaled up or down by 25%. Similarly to before, the barycenter of the triangle should not move due to the scaling. For this task, you can directly edit the position of the vertices before calling the rasterization function, or you can add a transformation in the vertex shader.

Ex.3: Colors
------------

Add the possibility to paint the color of each vertex in the scene. Color mode is enabled by the key <kbd>c</kbd>. In this mode, every mouse click will select the vertex closer to the current mouse position. After a vertex is selected, pressing a key from <kbd>1</kbd> to <kbd>9</kbd> will change its color (the colors that you use are not important, you can pick whatever colors you like). The color should be interpolated linearly inside the triangles.


Ex.4: Shader Translation/Scaling/Rotation
----------------------------------------------------

Instead of changing the coordinates of the vertices passed to the rasterizer, do all spatial transformations inside the vertex shader. While this will not make any difference in performance when using the software rasterizer, this will have a major impact when an hardware rasterizer is used.

Ex.5: View Control
------------------

Add to the application the capability of changing the camera. The following actions should be supported:

- <kbd>+</kbd> should increase the zoom by 20% zooming in in the center of the screen.

- <kbd>-</kbd> should decrease the zoom by 20% zooming out in the center of the screen.

- <kbd>w</kbd>, <kbd>a</kbd>, <kbd>s</kbd>, and <kbd>d</kbd> should pan the view by 20% of the visible part of the scene, i.e. translate the entire scene, respectively down, right, up and left by 20% of the window size.

This should NOT be implemented by changing the coordinates of the objects in the scene. You must add a view matrix to the vertex shader (as a uniform) that is transforming the position of the vertices of the triangles before they are rendered. Note that you will also have to transform the screen coordinates using the inverse of the view matrix, to ensure that the user interaction will adapt to the current view. ![image](img/view.png)


Ex.6: Add Keyframing
--------------------

### Linear Interpolation

Add the possibility to keyframe one property of an object (size, position, or rotation) and create an animation using linear interpolation between the keyframes. You can use a timer to make the animation automatic, or you could move to the next frame at the press of a button.

For example, you can press a key <kbd>K</kbd> to create a new frame of your animation. Move your triangles to the desired location, then hit <kbd>K</kbd> again to create a new frame. Once this is done, you can press <kbd>A</kbd> to play the animation. You can also add a shortcut <kbd>C</kbd> to clear the frames of the animation, keeping only the first frame. This is only a suggestion, but the essential to **not** hardcode the frame properties in your C++ code: we should be able to create our own animations on the spot.

### Bézier Curves

Add the option to use a [Bézier curve](https://en.wikipedia.org/wiki/B%C3%A9zier_curve) instead of a linear interpolation to animate the desired property (e.g. position) between each keyframe. For example, pressing <kbd>B</kbd> instead of <kbd>A</kbd> will play/advance the animation using the Bézier curve instead of a linear interpolation.
