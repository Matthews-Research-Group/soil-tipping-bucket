#ifndef SOIL_WATER_TILEDRAIN_H
#define SOIL_WATER_TILEDRAIN_H

#include "../framework/module.h"
#include "../framework/state_map.h"
#include "soil_water_flow_functions.h"

namespace BioCroWater 
{
/**
 * @class soil_water_downflow
 *
 * @brief Calculates water flow through the soil.
 *
 */
class soil_water_tiledrain : public direct_module
{
   public:
    soil_water_tiledrain(
        state_map const& input_quantities,
        state_map* output_quantities)
        : direct_module(),

          // Get references to input quantities
          tile_drain_depth{get_input(input_quantities, "tile_drain_depth")}, // FLDD. cm
          tile_drainage_rate{get_input(input_quantities, "tile_drainage_rate")}, // ETDR

          // Inputs for layer 1
          soil_depth_1{get_input(input_quantities, "soil_depth_1")},
          soil_water_content_1{get_input(input_quantities, "soil_water_content_1")},
          soil_field_capacity_1{get_input(input_quantities, "soil_field_capacity_1")},
          soil_saturation_capacity_1{get_input(input_quantities, "soil_saturation_capacity_1")},
          deltaS_1{get_input(input_quantities, "deltaS_1")},

          // Inputs for layer 2
          soil_depth_2{get_input(input_quantities, "soil_depth_2")},
          soil_water_content_2{get_input(input_quantities, "soil_water_content_2")},
          soil_field_capacity_2{get_input(input_quantities, "soil_field_capacity_2")},
          soil_saturation_capacity_2{get_input(input_quantities, "soil_saturation_capacity_2")},
          deltaS_2{get_input(input_quantities, "deltaS_2")},

          // Inputs for layer 3
          soil_depth_3{get_input(input_quantities, "soil_depth_3")},
          soil_water_content_3{get_input(input_quantities, "soil_water_content_3")},
          soil_field_capacity_3{get_input(input_quantities, "soil_field_capacity_3")},
          soil_saturation_capacity_3{get_input(input_quantities, "soil_saturation_capacity_3")},
          deltaS_3{get_input(input_quantities, "deltaS_3")},

          // Inputs for layer 4
          soil_depth_4{get_input(input_quantities, "soil_depth_4")},
          soil_water_content_4{get_input(input_quantities, "soil_water_content_4")},
          soil_field_capacity_4{get_input(input_quantities, "soil_field_capacity_4")},
          soil_saturation_capacity_4{get_input(input_quantities, "soil_saturation_capacity_4")},
          deltaS_4{get_input(input_quantities, "deltaS_4")},

          // Inputs for layer 5
          soil_depth_5{get_input(input_quantities, "soil_depth_5")},
          soil_water_content_5{get_input(input_quantities, "soil_water_content_5")},
          soil_field_capacity_5{get_input(input_quantities, "soil_field_capacity_5")},
          soil_saturation_capacity_5{get_input(input_quantities, "soil_saturation_capacity_5")},
          deltaS_5{get_input(input_quantities, "deltaS_5")},

          // Inputs for layer 6
          soil_depth_6{get_input(input_quantities, "soil_depth_6")},
          soil_water_content_6{get_input(input_quantities, "soil_water_content_6")},
          soil_field_capacity_6{get_input(input_quantities, "soil_field_capacity_6")},
          soil_saturation_capacity_6{get_input(input_quantities, "soil_saturation_capacity_6")},
          deltaS_6{get_input(input_quantities, "deltaS_6")},

          // Get pointers to output quantities
          td_layer_num_op{get_op(output_quantities, "td_layer_num")},
          cumulative_tile_flow_op{get_op(output_quantities, "cumulative_tile_flow")},
          total_tile_flow_op{get_op(output_quantities, "total_tile_flow")},

          head_op{get_op(output_quantities, "head")},
          tdf_avail_op{get_op(output_quantities, "tdf_avail")},
          topsat_op{get_op(output_quantities, "topsat")},
          tile_conductivity_op{get_op(output_quantities, "tile_conductivity")},

          deltaT_1_op{get_op(output_quantities, "deltaT_1")},
          deltaT_2_op{get_op(output_quantities, "deltaT_2")},
          deltaT_3_op{get_op(output_quantities, "deltaT_3")},
          deltaT_4_op{get_op(output_quantities, "deltaT_4")},
          deltaT_5_op{get_op(output_quantities, "deltaT_5")},
          deltaT_6_op{get_op(output_quantities, "deltaT_6")}

    {
    }
    static string_vector get_inputs();
    static string_vector get_outputs();
    static std::string get_name() { return "soil_water_tiledrain"; }

   private:
    // References to input quantities
    double const& tile_drain_depth;
    double const& tile_drainage_rate;

    // Inputs for layer 1
    double const& soil_depth_1;
    double const& soil_water_content_1;
    double const& soil_field_capacity_1;
    double const& soil_saturation_capacity_1;
    double const& deltaS_1;

    // Inputs for layer 2
    double const& soil_depth_2;
    double const& soil_water_content_2;
    double const& soil_field_capacity_2;
    double const& soil_saturation_capacity_2;
    double const& deltaS_2;

    // Inputs for layer 3
    double const& soil_depth_3;
    double const& soil_water_content_3;
    double const& soil_field_capacity_3;
    double const& soil_saturation_capacity_3;
    double const& deltaS_3;

    // Inputs for layer 4
    double const& soil_depth_4;
    double const& soil_water_content_4;
    double const& soil_field_capacity_4;
    double const& soil_saturation_capacity_4;
    double const& deltaS_4;

    // Inputs for layer 5
    double const& soil_depth_5;
    double const& soil_water_content_5;
    double const& soil_field_capacity_5;
    double const& soil_saturation_capacity_5;
    double const& deltaS_5;

    // Inputs for layer 6
    double const& soil_depth_6;
    double const& soil_water_content_6;
    double const& soil_field_capacity_6;
    double const& soil_saturation_capacity_6;
    double const& deltaS_6;

    // Pointers to output quantities
    double* td_layer_num_op;
    double* cumulative_tile_flow_op;
    double* total_tile_flow_op;

    double* head_op;
    double* tdf_avail_op;
    double* topsat_op;
    double* tile_conductivity_op;

    double* deltaT_1_op;
    double* deltaT_2_op;
    double* deltaT_3_op;
    double* deltaT_4_op;
    double* deltaT_5_op;
    double* deltaT_6_op;

    // Main operation
    void do_operation() const;
};

string_vector soil_water_tiledrain::get_inputs()
{
    return {
        "tile_drain_depth",   // Depth of the tile drain, (cm)
        "tile_drainage_rate", // Tile drainage rate 0.2 (from DSSAT)

        "soil_depth_1",
        "soil_water_content_1",
        "soil_field_capacity_1",
        "soil_saturation_capacity_1",
        "deltaS_1",         // Change in soil water content due to drainage in layer l (cm3 [water] / cm3 [soil])

        "soil_depth_2",
        "soil_water_content_2",
        "soil_field_capacity_2",
        "soil_saturation_capacity_2",
        "deltaS_2",

        "soil_depth_3",
        "soil_water_content_3",
        "soil_field_capacity_3",
        "soil_saturation_capacity_3",
        "deltaS_3",

        "soil_depth_4",
        "soil_water_content_4",
        "soil_field_capacity_4",
        "soil_saturation_capacity_4",
        "deltaS_4",

        "soil_depth_5",
        "soil_water_content_5",
        "soil_field_capacity_5",
        "soil_saturation_capacity_5",
        "deltaS_5",

        "soil_depth_6",
        "soil_water_content_6",
        "soil_field_capacity_6",
        "soil_saturation_capacity_6",
        "deltaS_6"
    };
}

string_vector soil_water_tiledrain::get_outputs()
{
    return {
        "td_layer_num",         // tile drain layer number
        "cumulative_tile_flow", // cm
        "total_tile_flow",
        "head",                 // cm
        "tdf_avail",            // soil water available to drain, cm 
        "topsat",               // to most saturated layer
        "tile_conductivity",    // cm/hr
        "deltaT_1",  // Change in soil water content due to tile drainage in layer 1 (cm3 [water] / cm3 [soil])
        "deltaT_2",  // Change in soil water content due to tile drainage in layer 2 (cm3 [water] / cm3 [soil])
        "deltaT_3",  // Change in soil water content due to tile drainage in layer 3 (cm3 [water] / cm3 [soil])
        "deltaT_4",  // Change in soil water content due to tile drainage in layer 4 (cm3 [water] / cm3 [soil])
        "deltaT_5",  // Change in soil water content due to tile drainage in layer 5 (cm3 [water] / cm3 [soil])
        "deltaT_6"   // Change in soil water content due to tile drainage in layer 6 (cm3 [water] / cm3 [soil])
    };
}

void soil_water_tiledrain::do_operation() const
{
    int nlayers = 6;

    double soil_depth[] = {
        soil_depth_1,
        soil_depth_2,
        soil_depth_3,
        soil_depth_4,
        soil_depth_5,
        soil_depth_6};

    double soil_water_content[] = {
        soil_water_content_1,
        soil_water_content_2,
        soil_water_content_3,
        soil_water_content_4,
        soil_water_content_5,
        soil_water_content_6};    

    double soil_field_capacity[] = {
        soil_field_capacity_1,
        soil_field_capacity_2,
        soil_field_capacity_3,
        soil_field_capacity_4,
        soil_field_capacity_5,
        soil_field_capacity_6};

    double soil_saturation_capacity[] = {
        soil_saturation_capacity_1,
        soil_saturation_capacity_2,
        soil_saturation_capacity_3,
        soil_saturation_capacity_4,
        soil_saturation_capacity_5,
        soil_saturation_capacity_6};

    double sw_delta_S[] = {
        deltaS_1,
        deltaS_2,
        deltaS_3,
        deltaS_4,
        deltaS_5,
        deltaS_6};


    int td_layer_num = 0;  // Layer number containing the tile drain
    if (tile_drain_depth <= 0.0) // Missing data or no tile
        td_layer_num = -99; 

    // Find layer number for tile
    else {
        double cumulative_depth = soil_depth[0];
        for (int l = 1; l < nlayers; l++){
            cumulative_depth = cumulative_depth + soil_depth[l];

            if ((cumulative_depth >= tile_drain_depth) &&  
                   ((cumulative_depth - soil_depth[l]) < tile_drain_depth)){
              
                td_layer_num = l; // Layer number containing the tile drain
            }
        }
    }
    tileDrain_str tileDrain;
    
    if (td_layer_num > 0){

        tileDrain = tile_flow(
            nlayers, 
            td_layer_num,
            tile_drainage_rate,
            soil_depth,
            soil_water_content,
            soil_field_capacity,
            soil_saturation_capacity, 
            sw_delta_S);
    // Skip tile drain section if depth of tile is <= 0
    } else {

        tileDrain.total_tile_flow = 0.0;
        tileDrain.cumulative_tile_flow = 0.0;
        for (int l = 0; l < nlayers; l++){
              tileDrain.sw_delta_T[l] = 0.0;
        }
    }
    // Update the output quantity list
    update(td_layer_num_op, td_layer_num);
    update(cumulative_tile_flow_op, tileDrain.cumulative_tile_flow);
    update(total_tile_flow_op, tileDrain.total_tile_flow);
    
    update(head_op, tileDrain.head);

    update(tdf_avail_op, tileDrain.tdf_avail);
    update(topsat_op, tileDrain.topsat);
    update(tile_conductivity_op, tileDrain.tile_drain_conductivity);

    update(deltaT_1_op, tileDrain.sw_delta_T[0]);
    update(deltaT_2_op, tileDrain.sw_delta_T[1]);
    update(deltaT_3_op, tileDrain.sw_delta_T[2]);
    update(deltaT_4_op, tileDrain.sw_delta_T[3]);
    update(deltaT_5_op, tileDrain.sw_delta_T[4]);
    update(deltaT_6_op, tileDrain.sw_delta_T[5]);
}

}  // namespace BioCroWater 
#endif
