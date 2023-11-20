import { NgModule } from '@angular/core';
import { BrowserModule } from '@angular/platform-browser';
import { AppComponent } from './app.component';
import { BrowserAnimationsModule } from '@angular/platform-browser/animations';
import { HttpClientModule } from '@angular/common/http';
import { BlockLibraryComponent } from './block-library/block-library.component';
import { CanvasComponent } from './canvas/canvas.component';
import { OutputDisplayComponent } from './output-display/output-display.component';
import { MatCardModule } from '@angular/material/card';
import { CodeBlockComponent } from './code-block/code-block.component';

@NgModule({
  declarations: [
    AppComponent,
    BlockLibraryComponent,
    CanvasComponent,
    OutputDisplayComponent,
    CodeBlockComponent,
  ],
  imports: [
    BrowserModule,
    BrowserAnimationsModule,
    HttpClientModule,
    MatCardModule,
  ],
  providers: [],
  bootstrap: [AppComponent]
})

export class AppModule { }
