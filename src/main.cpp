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
    GuiSetStyle(DEFAULT, TEXT_SIZE, 17);

    // Conversion factor from meters to pixels
    const float meter = std::min(screenWidth, screenHeight) * 0.1f;

    // Cart dimensions in meters
    const Vector2 cartDims = {meter, meter * 0.5f};
    const float wheelRadius = meter * 0.2f;
    // Center of cart body in pixels
    const Vector2 cartTrackCenter = {screenWidth * 0.75f / 2.0f - cartDims.x / 2.0f,
                                     screenHeight * 0.55f - cartDims.y - wheelRadius};

    // X, Y coordinates of DIP system in pixels
    Vector2 cartPos;
    Vector2 leftWheelPos;
    Vector2 rightWheelPos;
    Vector2 massOnePos;
    Vector2 massTwoPos;

    // Sliders (label, initial value, min value, max value)
    GUISlider sliders[] = {{"Cart Mass: %.2f kgs", 3.0f, 1.0f, 10.0f},
                           {"Mass 1: %.2f kgs", 1.0f, 0.5f, 5.0f},
                           {"Mass 2: %.2f kgs", 1.0f, 0.5f, 5.0f},
                           {"Length 1: %.2f m", 1.0f, 0.1f, 2.5f},
                           {"Length 2: %.2f m", 1.0f, 0.1f, 2.5f},
                           {"Cart Damping: %.2f Ns/m", 0.0f, 0.0f, 5.0f},
                           {"Theta 1 Damping: %.2f Nms/rad", 0.0f, 0.0f, 3.0f},
                           {"Theta 2 Damping: %.2f Nms/rad", 0.0f, 0.0f, 3.0f},
                           {"Theta 1: %.2f rad", 0.0f, -M_PIf, M_PIf},
                           {"Theta 2: %.2f rad", 0.0f, -M_PIf, M_PIf}};

    // Start button value
    bool startBtn = false;

    while (!WindowShouldClose()) {

        // Update postion of the DIP system
        cartPos = cartTrackCenter;  // While not started, cart stays at center
        leftWheelPos = {cartPos.x, cartPos.y + cartDims.y};
        rightWheelPos = {cartPos.x + cartDims.x, cartPos.y + cartDims.y};
        massOnePos = {cartPos.x + meter * sliders[3].currValue * sinf(sliders[8].currValue) + cartDims.x / 2.0f,
                      cartPos.y - meter * sliders[3].currValue * cosf(sliders[8].currValue) + cartDims.y / 2.0f};
        massTwoPos = {massOnePos.x + meter * sliders[4].currValue * sinf(sliders[8].currValue + sliders[9].currValue),
                      massOnePos.y - meter * sliders[4].currValue * cosf(sliders[8].currValue +  sliders[9].currValue)};

        BeginDrawing();

        ClearBackground(WHITE);

        // Cart Track
        DrawLineEx((Vector2){screenWidth * 0.05, screenHeight * 0.55},
                   (Vector2){screenWidth * 0.7, screenHeight * 0.55}, 
                   5,
                   Fade(BLACK, 0.8f));
        // Cart (Darker with increasing mass)
        DrawRectangleV(cartPos, cartDims,
                       Fade(BLUE, sliders[0].currValue / 20.0f + 0.5f));
        // Left wheel
        DrawCircleV(leftWheelPos, wheelRadius, Fade(RED, 0.95f));
        // Right wheel
        DrawCircleV(rightWheelPos, wheelRadius, Fade(RED, 0.95f));
        // Link 1
        DrawLineEx((Vector2){cartPos.x + cartDims.x / 2, cartPos.y + cartDims.y / 2},
                    massOnePos, 10, Fade(BLACK, 0.9f));
        // Link 2
        DrawLineEx(massOnePos, massTwoPos, 10, Fade(BLACK, 0.9f));
        // Mass 1 (Radius grows with mass)
        DrawCircleV(massOnePos, (sqrt(sliders[1].currValue) + 1.0) * meter / 10, RED);
        // Mass 2 (Radius grows with mass)
        DrawCircleV(massTwoPos, (sqrt(sliders[2].currValue) + 1.0) * meter / 10, RED);

        // Make GUI area on right quarter of window
        DrawLine(screenWidth * 0.75, 0, screenWidth * 0.75, screenWidth,
                 Fade(LIGHTGRAY, 0.6f));
        DrawRectangle(screenWidth * 0.75, 0, screenWidth * 0.25, screenHeight,
                      Fade(LIGHTGRAY, 0.3f));

        // Add GUI sliders
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
        if (GuiButton((Rectangle){screenWidth * 0.775, screenHeight * 0.675,
                                  screenWidth * 0.2, screenHeight * 0.025},
                       "Start")) {
            startBtn = true;
        }

        EndDrawing();
    }

    CloseWindow();

    return 0;
}