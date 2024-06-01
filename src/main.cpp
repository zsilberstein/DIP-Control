/**
 * @file main.cpp
 * @author Zach Silberstein (zach.silberstein@gmail.com)
 * @brief Dynamical simulation and control of the unstable double inverted
 * pendulum (DIP) on a cart.
 * @date 2024-05-31
 *
 * @copyright Copyright (c) 2024
 *
 */
#include <iostream>
#include "raylib.h"
#define RAYGUI_IMPLEMENTATION
#include "raygui.h"

/**
 * 
 * @struct GUISlider
 * @brief Stores data for a GUI slider.
 * @param label Label above the slider.
 * @param currVale The current value of the slider.
 * @param minValue The minumum value of the slider.
 * @param maxValue The maximum value of the slider.
 * 
 */
struct GUISlider {
    const char *label;
    float currValue;
    const float minValue;
    const float maxValue;
};

int main(void) {

    // Set up window
    const int screenWidth = 1600;
    const int screenHeight = 900;
    InitWindow(screenWidth, screenHeight, "DIP-Control");

    // Set style
    GuiSetStyle(DEFAULT, TEXT_SIZE, 17 * screenWidth / 1600);

    // sliders (label, initial value, min value, max value)
    GUISlider sliders[] = {{"Cart Mass: %.2f kgs", 3.0f, 1.0f, 10.0f},
                           {"Mass 1: %.2f kgs", 1.0f, 0.5f, 5.0f},
                           {"Mass 2: %.2f kgs", 1.0f, 0.5f, 5.0f},
                           {"Length 1: %.2f m", 1.0f, 0.1f, 2.5f},
                           {"Length 2: %.2f m", 1.0f, 0.1f, 2.5f},
                           {"Cart Damping: %.2f Ns/m", 0.0f, 0.0f, 5.0f},
                           {"Theta 1 Damping: %.2f Nms/rad", 0.0f, 0.0f, 3.0f},
                           {"Theta 2 Damping: %.2f Nms/rad", 0.0f, 0.0f, 3.0f}};

    // Start button value
    bool start = false;

    while (!WindowShouldClose()) {
        BeginDrawing();

        ClearBackground(WHITE);

        // Make GUI area on right quarter of window
        DrawLine(screenWidth * 0.75, 0, screenWidth * 0.75, screenWidth,
                 Fade(LIGHTGRAY, 0.6f));
        DrawRectangle(screenWidth * 0.75, 0, screenWidth * 0.25, screenHeight,
                      Fade(LIGHTGRAY, 0.3f));

        // Add sliders
        for (size_t i = 0; i < std::size(sliders); ++i) {
            GuiLabel((Rectangle){screenWidth * 0.775,
                                 screenHeight * (0.07f + i * 0.06f), 1400, 24},
                                 TextFormat(sliders[i].label, sliders[i].currValue));
            GuiSliderBar((Rectangle){screenWidth * 0.775,
                                     screenHeight * (0.09f + i * 0.06f),
                                     screenWidth * 0.2, screenHeight * 0.025},
                         NULL, NULL, &sliders[i].currValue, sliders[i].minValue,
                         sliders[i].maxValue);
        }

        // Add start button
        if (GuiButton((Rectangle){screenWidth * 0.775, screenHeight * 0.55,
                                  screenWidth * 0.2, screenHeight * 0.025},
                       "Start")) {
            start = true;
        }

        EndDrawing();
    }

    CloseWindow();

    return 0;
}