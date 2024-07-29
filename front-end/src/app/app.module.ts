import { AppComponent } from './app.component';
import { BrowserAnimationsModule } from '@angular/platform-browser/animations';
import { BrowserModule } from '@angular/platform-browser';
import { FormsModule } from '@angular/forms';
import { provideHttpClient, withInterceptorsFromDi } from '@angular/common/http';
import { MatButtonModule } from '@angular/material/button';
import { MatCardModule } from '@angular/material/card';
import { MatFormFieldModule } from '@angular/material/form-field';
import { MatProgressSpinnerModule } from '@angular/material/progress-spinner';
import { MatSelectModule } from '@angular/material/select';
import { MatSnackBarModule } from '@angular/material/snack-bar';
import { MatTooltipModule } from '@angular/material/tooltip';
import { NgModule } from '@angular/core';

import { BlockLibraryComponent } from './block-library/block-library.component';
import { CanvasComponent } from './canvas/canvas.component';
import { CodeBlockComponent } from './code-block/code-block.component';
import { OutputDisplayComponent } from './output-display/output-display.component';
import { SOCKET, socket } from './socket';

@NgModule({ declarations: [
  AppComponent,
  BlockLibraryComponent,
  CanvasComponent,
  CodeBlockComponent,
  OutputDisplayComponent,
],
bootstrap: [AppComponent], imports: [BrowserModule,
  BrowserAnimationsModule,
  FormsModule,
  MatButtonModule,
  MatCardModule,
  MatFormFieldModule,
  MatProgressSpinnerModule,
  MatSelectModule,
  MatSnackBarModule,
  MatTooltipModule], providers: [
  { provide: SOCKET, useValue: socket },
  provideHttpClient(withInterceptorsFromDi()),
] })

export class AppModule { }
