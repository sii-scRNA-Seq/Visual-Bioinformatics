
import { BrowserAnimationsModule } from '@angular/platform-browser/animations';
import { By } from '@angular/platform-browser';
import { ComponentFixture, TestBed } from '@angular/core/testing';
import { FormsModule } from '@angular/forms';
import { MatCardModule } from '@angular/material/card';
import { MatFormFieldModule } from '@angular/material/form-field';
import { MatIconModule } from '@angular/material/icon';
import { MatSelectModule } from '@angular/material/select';
import { MatSnackBarModule } from '@angular/material/snack-bar';

import { BlockService } from '../block.service';
import { CodeBlockComponent } from './code-block.component';
import { CurrentDatasetService } from '../current-dataset.service';
import { MockBlockService } from '../mock-block.service';
import { MockCurrentDatasetService } from '../mock-current-dataset.service';
import { MockOutputService } from '../mock-output.service';
import { OutputService } from '../output.service';

describe('CodeBlockComponent', () => {
  let component: CodeBlockComponent;
  let fixture: ComponentFixture<CodeBlockComponent>;

  beforeEach(() => {
    TestBed.configureTestingModule({
      declarations: [
        CodeBlockComponent
      ],
      imports: [
        BrowserAnimationsModule,
        FormsModule,
        MatCardModule,
        MatFormFieldModule,
        MatIconModule,
        MatSelectModule,
        MatSnackBarModule
      ],
      providers: [
        { provide: BlockService, useClass: MockBlockService },
        { provide: CurrentDatasetService, useClass: MockCurrentDatasetService },
        { provide: OutputService, useClass: MockOutputService }
      ]
    });
    fixture = TestBed.createComponent(CodeBlockComponent);
    component = fixture.componentInstance;
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });

  describe('Removing Blocks', () => {
    it('blockService.removeBlock should be called with loaddata when remove button is clicked', () => {
      component.block = {
        blockId: 'loaddata',
        blockUUID: '1',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [],
      };
      fixture.detectChanges();
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'removeBlock');
      const button = fixture.debugElement.query(By.css('button'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.removeBlock).toHaveBeenCalledOnceWith('1');
    });

    it('blockService.removeBlock should be called with basicfiltering when remove button is clicked', () => {
      component.block = {
        blockId: 'basicfiltering',
        blockUUID: '2',
        title: 'Basic Filtering',
        possibleChildBlocks: [],
        parameters: [],
      };
      fixture.detectChanges();
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'removeBlock');
      const button = fixture.debugElement.query(By.css('button'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.removeBlock).toHaveBeenCalledOnceWith('2');
    });

    it('blockService.removeBlock should be called with qcplots when remove button is clicked', () => {
      component.block = {
        blockId: 'qcplots',
        blockUUID: '3',
        title: 'Quality Control Plots',
        possibleChildBlocks: [],
        parameters: [],
      };
      fixture.detectChanges();
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'removeBlock');
      const button = fixture.debugElement.query(By.css('button'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.removeBlock).toHaveBeenCalledOnceWith('3');
    });

    it('blockService.removeBlock should be called with qcfiltering when remove button is clicked', () => {
      component.block = {
        blockId: 'qcfiltering',
        blockUUID: '4',
        title: 'Quality Control Filtering',
        possibleChildBlocks: [],
        parameters: [],
      };
      fixture.detectChanges();
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'removeBlock');
      const button = fixture.debugElement.query(By.css('button'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.removeBlock).toHaveBeenCalledOnceWith('4');
    });

    it('blockService.removeBlock should be called with variablegenes when remove button is clicked', () => {
      component.block = {
        blockId: 'variablegenes',
        blockUUID: '5',
        title: 'Identify Highly Variable Genes',
        possibleChildBlocks: [],
        parameters: [],
      };
      fixture.detectChanges();
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'removeBlock');
      const button = fixture.debugElement.query(By.css('button'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.removeBlock).toHaveBeenCalledOnceWith('5');
    });

    it('blockService.removeBlock should be called with pca when remove button is clicked', () => {
      component.block = {
        blockId: 'pca',
        blockUUID: '6',
        title: 'Principal Component Analysis',
        possibleChildBlocks: [],
        parameters: [],
      };
      fixture.detectChanges();
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'removeBlock');
      const button = fixture.debugElement.query(By.css('button'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.removeBlock).toHaveBeenCalledOnceWith('6');
    });

    it('blockService.removeBlock should be called with integration when remove button is clicked', () => {
      component.block = {
        blockId: 'integration',
        blockUUID: '7',
        title: 'Integration',
        possibleChildBlocks: [],
        parameters: [],
      };
      fixture.detectChanges();
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'removeBlock');
      const button = fixture.debugElement.query(By.css('button'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.removeBlock).toHaveBeenCalledOnceWith('7');
    });

    it('blockService.removeBlock should be called with runumap when remove button is clicked', () => {
      component.block = {
        blockId: 'runumap',
        blockUUID: '8',
        title: 'Run UMAP',
        possibleChildBlocks: [],
        parameters: [],
      };
      fixture.detectChanges();
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'removeBlock');
      const button = fixture.debugElement.query(By.css('button'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.removeBlock).toHaveBeenCalledOnceWith('8');
    });

    it('should be available when blocks are not being executed', () => {
      component.block = {
        blockId: 'loaddata',
        blockUUID: '',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [],
      };
      fixture.detectChanges();
      component.executingBlocks = false;
      fixture.detectChanges(); 
      expect(fixture.debugElement.query(By.css('button')).nativeElement.disabled).toEqual(false);
    });

    it('should become disabled while blocks are being executed', () => {
      component.block = {
        blockId: 'loaddata',
        blockUUID: '',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [],
      };
      fixture.detectChanges();
      component.executingBlocks = false;
      fixture.detectChanges(); 
      expect(fixture.debugElement.query(By.css('button')).nativeElement.disabled).toEqual(false);
      component.executingBlocks = true;
      fixture.detectChanges(); 
      expect(fixture.debugElement.query(By.css('button')).nativeElement.disabled).toEqual(true);
    });

    it('should become available once blocks have stopped being executed', () => {
      component.block = {
        blockId: 'loaddata',
        blockUUID: '',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [],
      };
      fixture.detectChanges();
      component.executingBlocks = true;
      fixture.detectChanges(); 
      expect(fixture.debugElement.query(By.css('button')).nativeElement.disabled).toEqual(true);
      component.executingBlocks = false;
      fixture.detectChanges(); 
      expect(fixture.debugElement.query(By.css('button')).nativeElement.disabled).toEqual(false);
    });
  });

  describe('Parameter Inputs', () => {
    it('should be available when blocks are not being executed', () => {
      component.block = {
        blockId: 'loaddata',
        blockUUID: '',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [
          {type: 'InputParameter', key: 'test_param', text: 'Test Parameter', value: 0},
        ],
      };
      fixture.detectChanges();
      component.executingBlocks = false;
      fixture.detectChanges();
      expect(fixture.debugElement.query(By.css('input')).nativeElement.disabled).toBe(false);
    });

    it('should become disabled while blocks are being executed', () => {
      component.block = {
        blockId: 'loaddata',
        blockUUID: '',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [
          {type: 'InputParameter', key: 'test_param', text: 'Test Parameter', value: 0},
        ],
      };
      fixture.detectChanges();
      component.executingBlocks = false;
      fixture.detectChanges();
      expect(fixture.debugElement.query(By.css('input')).nativeElement.disabled).toBe(false);
      component.executingBlocks = true;
      fixture.detectChanges(); 
      expect(fixture.debugElement.query(By.css('input')).nativeElement.disabled).toBe(true);
    });

    it('should become available once blocks have stopped being executed', () => {
      component.block = {
        blockId: 'loaddata',
        blockUUID: '',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [
          {type: 'InputParameter', key: 'test_param', text: 'Test Parameter', value: 0},
        ],
      };
      fixture.detectChanges();
      component.executingBlocks = true;
      fixture.detectChanges(); 
      expect(fixture.debugElement.query(By.css('input')).nativeElement.disabled).toBe(true);
      component.executingBlocks = false;
      fixture.detectChanges(); 
      expect(fixture.debugElement.query(By.css('input')).nativeElement.disabled).toBe(false);
    });
  });

  describe('Parameter Selects', () => {
    it('should be available when blocks are not being executed', () => {
      component.block = {
        blockId: 'loaddata',
        blockUUID: '',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [
          {type: 'SelectParameter', key: 'test_param', text: 'Test Parameter', value: 'option1', options: [
            {key: 'option1', text: 'Option 1'},
            {key: 'option2', text: 'Option 2'}
          ]},
        ],
      };
      fixture.detectChanges();
      component.executingBlocks = false;
      fixture.detectChanges();
      expect(fixture.debugElement.query(By.css('mat-select')).nativeElement.attributes.getNamedItem('ng-reflect-disabled').value).toEqual('false');
    });

    it('should become disabled while blocks are being executed', () => {
      component.block = {
        blockId: 'loaddata',
        blockUUID: '',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [
          {type: 'SelectParameter', key: 'test_param', text: 'Test Parameter', value: 'option1', options: [
            {key: 'option1', text: 'Option 1'},
            {key: 'option2', text: 'Option 2'}
          ]},
        ],
      };
      fixture.detectChanges();
      component.executingBlocks = false;
      fixture.detectChanges();
      expect(fixture.debugElement.query(By.css('mat-select')).nativeElement.attributes.getNamedItem('ng-reflect-disabled').value).toEqual('false');
      component.executingBlocks = true;
      fixture.detectChanges(); 
      expect(fixture.debugElement.query(By.css('mat-select')).nativeElement.attributes.getNamedItem('ng-reflect-disabled').value).toEqual('true');
    });

    it('should become available once blocks have stopped being executed', () => {
      component.block = {
        blockId: 'loaddata',
        blockUUID: '',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [
          {type: 'SelectParameter', key: 'test_param', text: 'Test Parameter', value: 'option1', options: [
            {key: 'option1', text: 'Option 1'},
            {key: 'option2', text: 'Option 2'}
          ]},
        ],
      };
      fixture.detectChanges();
      component.executingBlocks = true;
      fixture.detectChanges(); 
      expect(fixture.debugElement.query(By.css('mat-select')).nativeElement.attributes.getNamedItem('ng-reflect-disabled').value).toEqual('true');
      component.executingBlocks = false;
      fixture.detectChanges(); 
      expect(fixture.debugElement.query(By.css('mat-select')).nativeElement.attributes.getNamedItem('ng-reflect-disabled').value).toEqual('false');
    });

    it('should call currentDatasetService.setCurrentDataset() when value changes and blockId is loaddata', () => {
      component.block = {
        blockId: 'loaddata',
        blockUUID: '',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [
          {type: 'SelectParameter', key: 'test_param', text: 'Test Parameter', value: 'option1', options: [
            {key: 'option1', text: 'Option 1'},
            {key: 'option2', text: 'Option 2'}
          ]},
        ],
      };
      fixture.detectChanges();
      const event = {value: 'option2'};
      const currentDatasetService: CurrentDatasetService = TestBed.inject(CurrentDatasetService);
      spyOn(currentDatasetService, 'setCurrentDataset');
      fixture.debugElement.query(By.css('mat-select')).triggerEventHandler('selectionChange', event);
      fixture.detectChanges();
      expect(currentDatasetService.setCurrentDataset).toHaveBeenCalledTimes(1);
    });

    it('should not call currentDatasetService.setCurrentDataset() when value changes and blockId is not loaddata', () => {
      component.block = {
        blockId: 'qcfiltering',
        blockUUID: '',
        title: 'Quality Control Filtering',
        possibleChildBlocks: [],
        parameters: [
          {type: 'SelectParameter', key: 'test_param', text: 'Test Parameter', value: 'option1', options: [
            {key: 'option1', text: 'Option 1'},
            {key: 'option2', text: 'Option 2'}
          ]},
        ],
      };
      fixture.detectChanges();
      const event = {value: 'option2'};
      const currentDatasetService: CurrentDatasetService = TestBed.inject(CurrentDatasetService);
      spyOn(currentDatasetService, 'setCurrentDataset');
      fixture.debugElement.query(By.css('mat-select')).triggerEventHandler('selectionChange', event);
      fixture.detectChanges();
      expect(currentDatasetService.setCurrentDataset).toHaveBeenCalledTimes(0);
    });
  });
});
