﻿#include "glfw-3.2.1/include/GLFW/glfw3.h"
#include "../SHtest/SH.h"

#include <iostream>
#include <string>
#include <cstdio>

#pragma warning(disable: 4996)

sh::SHS main_sh = std::vector<float>{1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0, 0, 0, 0};
float sh_pos_z = -2;
float sh_rot_x = -60;
float sh_rot_z = -120;
int cur_pos = 0;

GLFWwindow *main_win = nullptr;
int win_width = 512;
int win_height = 512;
bool mouse_button_l = false;
bool mouse_button_r = false;
int mouse_last_x = 0;
int mouse_last_y = 0;

void draw_sh(const sh::SHS &s, std::function<void(float, float, float)> colorFunc = [](float,float,float){});

int main(int argc, const char *argv[]) {
	glfwInit();
	main_win = glfwCreateWindow(win_width, win_height, "SHtest3D", nullptr, nullptr);

	glfwSetKeyCallback(
		main_win,
		[](GLFWwindow *window, int key, int scancode, int action, int mods){
			if (action != GLFW_PRESS) return;

			if (key == GLFW_KEY_ESCAPE) glfwSetWindowShouldClose(main_win, GL_TRUE);

			const bool key_left = key == GLFW_KEY_LEFT;
			const bool key_right = key == GLFW_KEY_RIGHT;
			const bool key_up = key == GLFW_KEY_UP;
			const bool key_down = key == GLFW_KEY_DOWN;
			const bool key_shift = mods & GLFW_MOD_SHIFT;

			if (key_shift) {
				if (key_left) main_sh = rotate(main_sh, -0.1, 0, 0);
				if (key_right) main_sh = rotate(main_sh, 0.1, 0, 0);
				if (key_up) main_sh = rotate(main_sh, 0, -0.1, 0);
				if (key_down) main_sh = rotate(main_sh, 0, 0.1, 0);
			} else {
				if (key_left && cur_pos > 0) --cur_pos;
				if (key_right) {
					++cur_pos;
					if (cur_pos >= main_sh.num_coeffs()) {
						sh::SHS new_sh(main_sh.max_l() + 1);
						new_sh.set_coeffs(main_sh.coeffs(), main_sh.num_coeffs());
						main_sh = new_sh;
					}
				}
				if (key_up || key_down) {
					std::vector<float> coeffs(main_sh.coeffs(), main_sh.coeffs() + main_sh.num_coeffs());
					coeffs[cur_pos] += key_up ? 0.1 : -0.1;
					main_sh.set_coeffs(coeffs);
				}
			}

			std::string title = "SHtest3D ";
			for (int c = 0; c < main_sh.num_coeffs(); ++c) {
				char buf[16];
				std::sprintf(buf, "%02.1f", main_sh.coeffs()[c]);
				if (c > 0) title += ", ";
				if (c == cur_pos) title += '[';
				title += buf;
				if (c == cur_pos) title += ']';
			}
			glfwSetWindowTitle(main_win, title.c_str());
		}
	);

	glfwSetWindowSizeCallback(
		main_win,
		[](GLFWwindow *window, int width, int height) {
			win_width = width;
			win_height = height;
		}
	);

	glfwSetMouseButtonCallback(
		main_win,
		[](GLFWwindow *window, int button, int action, int mods) {
			if (button == GLFW_MOUSE_BUTTON_LEFT) mouse_button_l = action == GLFW_PRESS;
			if (button == GLFW_MOUSE_BUTTON_RIGHT) mouse_button_r = action == GLFW_PRESS;
		}
	);

	glfwSetCursorPosCallback(
		main_win,
		[](GLFWwindow *window, double xpos, double ypos) {
			if (mouse_button_l) {
				sh_rot_z += xpos - mouse_last_x;
				sh_rot_x += ypos - mouse_last_y;
			}
			if (mouse_button_r) {
				sh_pos_z += (ypos - mouse_last_y) * 0.01f;
			}
			mouse_last_x = xpos;
			mouse_last_y = ypos;
		}
	);

	glfwMakeContextCurrent(main_win);

	while (!glfwWindowShouldClose(main_win)) {
		glfwWaitEvents();

		glViewport(0, 0, win_width, win_height);
		glClearColor(0.1f, 0.2f, 0.3f, 0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glEnable(GL_DEPTH_TEST);

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		const float a = win_width < win_height ? win_width : win_height;
		if (a <= 0) continue;
		const float n = 0.001f;
		const float w = win_width / a * n;
		const float h = win_height / a * n;
		glFrustum(-w, w, -h, h, n, 10);

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glTranslatef(0, 0, sh_pos_z);
		glRotatef(sh_rot_x, 1, 0, 0);
		glRotatef(sh_rot_z, 0, 0, 1);


		glColor3f(0.25f, 0.5f, 0.5f);
		draw_sh(main_sh, glColor3f);

		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glPolygonOffset(-0.1, -1);
		glEnable(GL_POLYGON_OFFSET_LINE);
		glColor3f(1, 1, 1);
		draw_sh(main_sh);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		glBegin(GL_LINES);
		glColor3f(1, 0, 0); glVertex3f(0, 0, 0); glVertex3f(1, 0, 0);
		glColor3f(0, 1, 0); glVertex3f(0, 0, 0); glVertex3f(0, 1, 0);
		glColor3f(0, 0, 1); glVertex3f(0, 0, 0); glVertex3f(0, 0, 1);
		glEnd();

		GLenum err = glGetError();
		if (err != GL_NO_ERROR) std::cout << "GL Error: " << err << std::endl;

		glfwSwapBuffers(main_win);
	}
	glfwDestroyWindow(main_win);
	glfwTerminate();
}

void draw_sh(const sh::SHS &s, std::function<void(float,float,float)> colorFunc) {
	glBegin(GL_QUADS);
	const int nDivT = 32;
	const int nDivP = 64;
	const float dt = 3.141592f / nDivT;
	const float dp = 2 * 3.141592f / nDivP;
	for (int t = 0; t < nDivT; ++t) {
		for (int p = 0; p < nDivP; ++p) {
			const float t0 = dt * t;
			const float p0 = dp * p;
			const float t1 = t0 + dt;
			const float p1 = p0 + dp;
			const float r0 = s.eval(t0, p0);
			const float r1 = s.eval(t1, p0);
			const float r2 = s.eval(t1, p1);
			const float r3 = s.eval(t0, p1);
			using namespace std;
			if (r0 > 0) colorFunc(1,0,0); else colorFunc(0,1,0);
			glVertex3f(r0*r0 * sin(t0) * cos(p0), r0*r0 * std::sin(t0) * sin(p0), r0*r0 * cos(t0));
			if (r1 > 0) colorFunc(1, 0, 0); else colorFunc(0, 1, 0);
			glVertex3f(r1*r1 * sin(t1) * cos(p0), r1*r1 * std::sin(t1) * sin(p0), r1*r1 * cos(t1));
			if (r2 > 0) colorFunc(1, 0, 0); else colorFunc(0, 1, 0);
			glVertex3f(r2*r2 * sin(t1) * cos(p1), r2*r2 * std::sin(t1) * sin(p1), r2*r2 * cos(t1));
			if (r3 > 0) colorFunc(1, 0, 0); else colorFunc(0, 1, 0);
			glVertex3f(r3*r3 * sin(t0) * cos(p1), r3*r3 * std::sin(t0) * sin(p1), r3*r3 * cos(t0));

		}
	}
	glEnd();
}
