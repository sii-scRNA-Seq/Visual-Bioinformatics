import { ComponentFixture, TestBed } from '@angular/core/testing';
import { DomSanitizer } from '@angular/platform-browser';
import { HttpClientTestingModule } from '@angular/common/http/testing';
import { MatCardModule } from '@angular/material/card';
import { MatSnackBarModule } from '@angular/material/snack-bar';

import { Block } from '../block.interface';
import { MockOutputService } from '../mock-output.service';
import { OutputDisplayComponent } from './output-display.component';
import { OutputService } from '../output.service';

describe('OutputDisplayComponent', () => {
  let component: OutputDisplayComponent;
  let fixture: ComponentFixture<OutputDisplayComponent>;


  beforeEach(() => {
    TestBed.configureTestingModule({
      declarations: [OutputDisplayComponent],
      imports: [
        HttpClientTestingModule,
        MatCardModule,
        MatSnackBarModule,
      ],
      providers: [
        { provide: OutputService, useClass: MockOutputService }
      ],
    });
    fixture = TestBed.createComponent(OutputDisplayComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });

  describe('OutputList', () => {
    it('should update when OutputService.outputs updates', () => {
      expect(component.outputList.length).toBe(0);
      const block: Block = {
        blockId: 'loaddata',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [],
      };
      const outputService: OutputService = TestBed.inject(OutputService);
      outputService.executeBlock(block);
      expect(component.outputList.length).toBe(1);
    });
  });
});
